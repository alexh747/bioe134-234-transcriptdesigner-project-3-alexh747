# This script is just to try out different transcript designers. The real improved transcript designer that
# I want to submit for project 3 is transcript_designer.py

# Code after the forbidden site and hairpin iteration modifications

import random
import csv
from collections import defaultdict, Counter

from genedesign.rbs_chooser import RBSChooser
from genedesign.models.transcript import Transcript
from genedesign.checkers.codon_checker import CodonChecker
from genedesign.checkers.forbidden_sequence_checker import ForbiddenSequenceChecker
from genedesign.checkers.internal_promoter_checker import PromoterChecker
from genedesign.checkers.hairpin_checker import hairpin_checker
from genedesign.checkers.gc_content_checker import GCContentChecker


class TranscriptDesigner:
    """
    Reverse translates a protein sequence into a codon-optimized DNA sequence
    and chooses an RBS. The design is iteratively refined using a stochastic 
    optimizer to satisfy all UCB_BioE134 quality checks.
    """

    def __init__(self):
        self.rbsChooser = None
        self.codonChecker = None
        self.forbiddenChecker = None
        self.promoterChecker = None
        self.gcChecker = None
        
        self.aa_to_codons = {}
        self.rare_codons = set()
        self._good_codons = {}

    def initiate(self) -> None:
        """Loads codon usage data and initializes every sub-component."""
        self.rbsChooser = RBSChooser()
        self.rbsChooser.initiate()

        self.codonChecker = CodonChecker()
        self.codonChecker.initiate()

        self.forbiddenChecker = ForbiddenSequenceChecker()
        self.forbiddenChecker.initiate()

        self.promoterChecker = PromoterChecker()
        self.promoterChecker.initiate()

        self.gcChecker = GCContentChecker()
        self.gcChecker.initiate()

        # Load codon usage data
        codon_usage_file = 'genedesign/data/codon_usage.txt'
        raw = defaultdict(list)
        try:
            with open(codon_usage_file, 'r') as f:
                reader = csv.reader(f, delimiter='\t')
                for row in reader:
                    if len(row) < 4: continue
                    codon = row[0].strip().upper()
                    aa = row[1].strip().upper()
                    freq = float(row[2].strip())
                    raw[aa].append((codon, freq))
        except FileNotFoundError:
            # Fallback if file isn't found during testing
            pass 

        self.aa_to_codons = {aa: sorted(clist, key=lambda x: -x[1]) for aa, clist in raw.items()}
        
        # Determine rare codons (frequency < 0.10)
        self.rare_codons = set()
        for aa, clist in self.aa_to_codons.items():
            for codon, freq in clist:
                if freq < 0.10:
                    self.rare_codons.add(codon)

        # Pre-compute allowed (non-rare) codons
        for aa, clist in self.aa_to_codons.items():
            good = [c for c, w in clist if c not in self.rare_codons]
            if not good:
                good = [c for c, w in clist] # Fallback if all are "rare"
            self._good_codons[aa] = good

        # Ensure stop codon exists
        if '*' not in self._good_codons:
            self._good_codons['*'] = ['TAA', 'TGA', 'TAG']

    def _get_random_good_codon(self, aa: str) -> str:
        """Picks a valid, non-rare codon uniformly to boost diversity."""
        options = self._good_codons.get(aa, ["NNN"])
        return random.choice(options)

    def _score_sequence(self, codons: list, rbs_utr: str) -> int:
        """
        Calculates a penalty score. 0 means it perfectly passes all checkers.
        """
        penalty = 0
        cds_str = ''.join(codons)
        transcript_dna = rbs_utr.upper() + cds_str

        # 1. Codon Usage Check (Needs diversity >= 0.5, rare <= 3, CAI >= 0.2)
        codons_pass, diversity, rare_count, cai = self.codonChecker.run(codons)
        if not codons_pass:
            if diversity < 0.5:
                penalty += int((0.5 - diversity) * 200) + 10
            if rare_count > 3:
                penalty += (rare_count - 3) * 20
            if cai < 0.2:
                penalty += 50
        
        # 2. GC Content Check
        gc_pass, gc_val, gc_msg = self.gcChecker.run(cds_str)
        if not gc_pass:
            penalty += 50

        # 3. Forbidden Sequence Check
        forb_pass, forb_site = self.forbiddenChecker.run(transcript_dna)
        if not forb_pass:
            penalty += 200

        # 4. Promoter Check
        prom_pass, prom_site = self.promoterChecker.run(transcript_dna)
        if not prom_pass:
            penalty += 150

        # 5. Hairpin Check
        hairpin_pass, hp_str = hairpin_checker(transcript_dna)
        if not hairpin_pass:
            penalty += 100

        return penalty

    def run(self, peptide: str, ignores: set) -> Transcript:
        """
        Main execution loop.
        """
        # 1. Formatting
        protein = peptide.upper().strip()
        if not protein.endswith('*'):
            protein += '*'
            
        # 2. Build initial diverse CDS
        codons = []
        for aa in protein:
            codons.append(self._get_random_good_codon(aa))
        
        # Get RBS
        cds_str = ''.join(codons)
        selectedRBS = self.rbsChooser.run(cds_str, ignores)
        rbs_utr = selectedRBS.utr

        # 3. Stochastic Hill Climbing Optimization Loop
        best_codons = list(codons)
        best_penalty = self._score_sequence(best_codons, rbs_utr)

        max_iterations = 1000
        
        for i in range(max_iterations):
            if best_penalty == 0:
                break # Perfect sequence found!

            test_codons = list(best_codons)

            # Mutate 1 to 3 random codons (excluding the very last stop codon if possible)
            num_mutations = random.randint(1, 3)
            for _ in range(num_mutations):
                idx = random.randint(0, len(test_codons) - 2)
                aa = protein[idx]
                test_codons[idx] = self._get_random_good_codon(aa)

            # Update RBS slightly just in case the first codons changed dramatically
            test_cds_str = ''.join(test_codons)
            
            test_penalty = self._score_sequence(test_codons, rbs_utr)

            # Accept if strictly better, or occasionally if equal (to escape flat local minima)
            if test_penalty < best_penalty or (test_penalty == best_penalty and random.random() < 0.2):
                best_codons = list(test_codons)
                best_penalty = test_penalty

        # Final packaging
        final_cds_str = ''.join(best_codons)
        # Re-run RBS chooser one last time to ensure it perfectly matches the 5' end of the CDS
        final_rbs = self.rbsChooser.run(final_cds_str, ignores)

        return Transcript(final_rbs, protein, best_codons)

if __name__ == "__main__":
    peptide = "MYPFIRTARMTVCAKKHVHLTRDAAEQLLADIDRRLDQLLPVEGERDVVGAAMREGALAPGKRIRPMLLLLTARDLGCAVSHDGLLDLACAVEMVHAASLILDDMPCMDDAKLRRGRPTIHSHYGEHVAILAAVALLSKAFGVIADADGLTPLAKNRAVSELSNAIGMQGLVQGQFKDLSEGDKPRSAEAILMTNHFKTSTLFCASMQMASIVANASSEARDCLHRFSLDLGQAFQLLDDLTDGMTDTGKDSNQDAGKSTLVNLLGPRAVEERLRQHLQLASEHLSAACQHGHATQHFIQAWFDKKLAAVS"

    designer = TranscriptDesigner()
    designer.initiate()
    ignores = set()
    transcript = designer.run(peptide, ignores)

    cds = ''.join(transcript.codons)
    print(f"CDS length: {len(cds)}")
    print(f"CDS ends with stop:  {cds[-3:] in ('TAA','TGA','TAG')}")
    print("Design complete.")
