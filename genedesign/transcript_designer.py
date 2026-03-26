# Code before the forbidden site and haiprin iteration modification

# New Version of Transcript Designer 
import random
import csv
from collections import defaultdict

from genedesign.rbs_chooser import RBSChooser
from genedesign.models.transcript import Transcript
from genedesign.checkers.codon_checker import CodonChecker
from genedesign.checkers.forbidden_sequence_checker import ForbiddenSequenceChecker
from genedesign.checkers.internal_promoter_checker import PromoterChecker
from genedesign.checkers.hairpin_checker import hairpin_checker
from genedesign.checkers.gc_content_checker import GCContentChecker
from genedesign.seq_utils.hairpin_counter import hairpin_counter
from genedesign.seq_utils.reverse_complement import reverse_complement


class TranscriptDesigner:
    """
    Reverse translates a protein sequence into a codon-optimized DNA sequence
    and chooses an RBS.  The design is iteratively refined until it passes all
    quality checks:

        * High Codon Adaptation Index (CAI) with adequate codon diversity
        * No forbidden restriction-enzyme sites or problematic motifs
        * No internal sigma-70 promoter sequences
        * Low hairpin count
        * Acceptable GC content (global and local)
    """

    def __init__(self):
        self.aa_to_codons = {}       # aa -> [(codon, freq), ...]
        self.codon_to_aa = {}        # codon -> aa
        self.rbsChooser = None
        self.codonChecker = None
        self.forbiddenChecker = None
        self.promoterChecker = None
        self.gcChecker = None
        # Pre-computed codon selection tables (populated in initiate)
        self._good_codons = {}       # aa -> [codon, ...]
        self._good_weights = {}      # aa -> [weight, ...]
        self._alt_codons = {}        # codon -> ([alt_codons], [weights])

    # ------------------------------------------------------------------
    # Initialisation
    # ------------------------------------------------------------------
    def initiate(self) -> None:
        """
        Loads codon usage data and initialises every sub-component
        (RBS chooser and all sequence checkers).
        """
        # --- RBS chooser ---
        self.rbsChooser = RBSChooser()
        self.rbsChooser.initiate()

        # --- Checkers ---
        self.codonChecker = CodonChecker()
        self.codonChecker.initiate()

        self.forbiddenChecker = ForbiddenSequenceChecker()
        self.forbiddenChecker.initiate()

        self.promoterChecker = PromoterChecker()
        self.promoterChecker.initiate()

        self.gcChecker = GCContentChecker()
        self.gcChecker.initiate()

        # --- Build amino-acid -> codon table from codon_usage.txt ---
        codon_usage_file = 'genedesign/data/codon_usage.txt'
        raw = defaultdict(list)
        with open(codon_usage_file, 'r') as f:
            reader = csv.reader(f, delimiter='\t')
            for row in reader:
                if len(row) < 4:
                    continue
                codon = row[0].strip()
                aa = row[1].strip()
                freq = float(row[2].strip())
                if aa == '*':
                    continue  # skip stop codons for normal selection
                raw[aa].append((codon, freq))
                self.codon_to_aa[codon] = aa

        # Store codons sorted by descending frequency for each amino acid
        self.aa_to_codons = {}
        for aa, codon_list in raw.items():
            self.aa_to_codons[aa] = sorted(codon_list, key=lambda x: -x[1])

        # Build a set of rare codons (freq < 0.10) to avoid during selection
        self.rare_codons = set()
        for aa, codon_list in self.aa_to_codons.items():
            for codon, freq in codon_list:
                if freq < 0.10:
                    self.rare_codons.add(codon)

        # Pre-compute filtered codon lists and weights for fast selection
        for aa, codon_list in self.aa_to_codons.items():
            good = [(c, w) for c, w in codon_list if c not in self.rare_codons]
            if not good:
                good = codon_list
            self._good_codons[aa] = [c for c, _ in good]
            self._good_weights[aa] = [w for _, w in good]

        # Pre-compute alternate codon tables (excluding each codon itself)
        for aa, codon_list in self.aa_to_codons.items():
            for codon, _ in codon_list:
                alts = [(c, w) for c, w in codon_list if c not in self.rare_codons and c != codon]
                if not alts:
                    alts = [(c, w) for c, w in codon_list if c != codon]
                if not alts:
                    self._alt_codons[codon] = ([codon], [1.0])
                else:
                    self._alt_codons[codon] = ([c for c, _ in alts], [w for _, w in alts])

    # ------------------------------------------------------------------
    # Codon selection helpers
    # ------------------------------------------------------------------
    def _weighted_random_codon(self, aa: str) -> str:
        """Pick a codon for *aa* with probability proportional to its usage
        frequency, excluding rare codons when possible."""
        return random.choices(self._good_codons[aa], self._good_weights[aa])[0]

    def _alternate_codon(self, codon: str) -> str:
        """Return a different synonymous codon chosen by weighted random selection,
        avoiding rare codons when possible.
        If the amino acid has only one codon (M, W) the same codon is returned."""
        codons_list, weights = self._alt_codons.get(codon, ([codon], [1.0]))
        return random.choices(codons_list, weights)[0]

    # ------------------------------------------------------------------
    # Initial CDS generation
    # ------------------------------------------------------------------
    def _initial_codons(self, peptide: str) -> list:
        """
        Generate an initial codon list using frequency-weighted random selection.
        The first codon is always ATG (start codon) regardless of the first amino acid.
        This naturally produces high CAI and good codon diversity because
        common codons are selected proportionally to their usage in E. coli.
        """
        codons = []
        for idx, aa in enumerate(peptide):
            if idx == 0:
                codons.append("ATG")  # always start with ATG
            else:
                codons.append(self._weighted_random_codon(aa))
        codons.append("TAA")  # stop codon
        return codons

    # ------------------------------------------------------------------
    # Checker wrappers
    # ------------------------------------------------------------------
    def _check_all(self, codons: list, rbs_utr: str) -> dict:
        """
        Run all checkers on the current codon list and return a dict of results.
        Each value is a tuple whose first element is a bool (True = passed).
        """
        cds = ''.join(codons)
        transcript_dna = rbs_utr.upper() + cds

        codon_result = self.codonChecker.run(codons)                  # (bool, div, rare, cai)
        forbidden_result = self.forbiddenChecker.run(transcript_dna)  # (bool, site)
        promoter_result = self.promoterChecker.run(transcript_dna)    # (bool, seq)
        hairpin_result = hairpin_checker(transcript_dna)              # (bool, str)
        gc_result = self.gcChecker.run(cds)                           # (bool, gc, msg)

        return {
            'codon':     codon_result,
            'forbidden': forbidden_result,
            'promoter':  promoter_result,
            'hairpin':   hairpin_result,
            'gc':        gc_result,
        }

    @staticmethod
    def _all_pass(results: dict) -> bool:
        """Return True if every checker passed."""
        return all(v[0] for v in results.values())

    @staticmethod
    def _score(results: dict) -> int:
        """Return the number of checkers that passed (higher is better)."""
        return sum(1 for v in results.values() if v[0])

    # ------------------------------------------------------------------
    # Targeted repair strategies
    # ------------------------------------------------------------------
    def _find_site_codon_range(self, codons: list, rbs_utr: str, site: str) -> tuple:
        """Locate *site* in the transcript DNA (forward or reverse complement)
        and return the codon index range (start, end) within the CDS that
        overlaps with it.  Returns None if not found on either strand."""
        cds = ''.join(codons)
        transcript_dna = rbs_utr.upper() + cds
        utr_len = len(rbs_utr)

        pos = transcript_dna.upper().find(site.upper())
        if pos == -1:
            # Check reverse complement strand
            rc = reverse_complement(transcript_dna)
            rc_pos = rc.upper().find(site.upper())
            if rc_pos == -1:
                return None
            pos = len(transcript_dna) - rc_pos - len(site)

        cds_start = max(0, pos - utr_len)
        cds_end = min(len(cds), cds_start + len(site) + 3)
        codon_start = max(0, cds_start // 3 - 1)
        codon_end = min(len(codons) - 1, cds_end // 3 + 1)
        return codon_start, codon_end

    def _fix_at_site(self, codons: list, rbs_utr: str, site: str) -> None:
        """Swap codons overlapping a problematic site to break it.
        Used for both forbidden sequences and internal promoters."""
        if site is None:
            return
        rng = self._find_site_codon_range(codons, rbs_utr, site)
        if rng is None:
            self._random_swap(codons, 3)
        else:
            for i in range(rng[0], rng[1]):
                codons[i] = self._alternate_codon(codons[i])

    def _random_swap(self, codons: list, n: int) -> list:
        """Swap *n* random coding positions (not the stop codon)."""
        n_coding = len(codons) - 1
        if n_coding <= 0:
            return codons
        n = min(n, n_coding)
        indices = random.sample(range(n_coding), n)
        for i in indices:
            codons[i] = self._alternate_codon(codons[i])
        return codons

    def _sweep_hairpins(self, codons: list, rbs_utr: str, max_passes: int = 10) -> list:
        """Walk through the transcript doing targeted codon swaps at every
        chunk that has more than 1 hairpin.  Repeat up to *max_passes* times."""
        chunk_size = 50
        overlap = 25

        for _pass in range(max_passes):
            cds = ''.join(codons)
            transcript_dna = rbs_utr.upper() + cds
            utr_len = len(rbs_utr)
            any_bad = False

            for i in range(0, len(transcript_dna) - chunk_size + 1, overlap):
                chunk = transcript_dna[i:i + chunk_size]
                hcount, _ = hairpin_counter(chunk, 3, 4, 9)
                if hcount > 1:
                    any_bad = True
                    # Map chunk region to codon indices
                    cds_start = max(0, i - utr_len)
                    cds_end = min(len(cds), i + chunk_size - utr_len)
                    codon_start = max(0, cds_start // 3)
                    codon_end = min(len(codons) - 1, cds_end // 3 + 1)

                    for ci in range(codon_start, codon_end):
                        codons[ci] = self._alternate_codon(codons[ci])

                    # Re-build transcript_dna after changes for next chunk
                    cds = ''.join(codons)
                    transcript_dna = rbs_utr.upper() + cds

            if not any_bad:
                break

        return codons

    # ------------------------------------------------------------------
    # Main entry point
    # ------------------------------------------------------------------
    def run(self, peptide: str, ignores: set) -> Transcript:
        """
        Designs a codon-optimized CDS for the given peptide and selects an RBS.

        Algorithm:
          1. Generate an initial CDS using frequency-weighted random codon selection.
          2. Pick an RBS.
          3. Single repair cycle: fix forbidden / promoters / hairpins.
          4. Track the best candidate across restarts if needed.
          5. Return the best (or first perfect) result.

        Parameters:
            peptide (str): The amino-acid sequence (single-letter codes, no stop).
            ignores (set): RBS options to skip.

        Returns:
            Transcript: The designed transcript with RBS and codon list.
        """
        MAX_RESTARTS = 2

        best_codons = None
        best_score = -1
        best_rbs = None

        for _restart in range(MAX_RESTARTS):
            # Step 1: initial codon selection
            codons = self._initial_codons(peptide)

            # Step 2: build CDS and choose RBS
            cds = ''.join(codons)
            selectedRBS = self.rbsChooser.run(cds, ignores)
            rbs_utr = selectedRBS.utr

            # Step 3: single repair cycle (no double checker calls)

            # Fix forbidden sites
            for _ in range(6):
                passed, site = self.forbiddenChecker.run(rbs_utr.upper() + ''.join(codons))
                if passed:
                    break
                self._fix_at_site(codons, rbs_utr, site)

            # Fix promoters
            for _ in range(6):
                passed, seq = self.promoterChecker.run(rbs_utr.upper() + ''.join(codons))
                if passed:
                    break
                self._fix_at_site(codons, rbs_utr, seq)

            # Hairpin sweep — only if needed
            passed_h, _ = hairpin_checker(rbs_utr.upper() + ''.join(codons))
            if not passed_h:
                codons = self._sweep_hairpins(codons, rbs_utr, max_passes=3)

            # Cleanup: re-fix forbidden/promoter that hairpin sweep may have introduced
            for _ in range(4):
                passed, site = self.forbiddenChecker.run(rbs_utr.upper() + ''.join(codons))
                if passed:
                    break
                self._fix_at_site(codons, rbs_utr, site)

            for _ in range(4):
                passed, seq = self.promoterChecker.run(rbs_utr.upper() + ''.join(codons))
                if passed:
                    break
                self._fix_at_site(codons, rbs_utr, seq)

            # Score this candidate
            results = self._check_all(codons, rbs_utr)
            if self._all_pass(results):
                return Transcript(selectedRBS, peptide, codons)

            score = self._score(results)
            if score > best_score:
                best_score = score
                best_codons = list(codons)
                best_rbs = selectedRBS

        # Return the best candidate found across all restarts
        return Transcript(best_rbs, peptide, best_codons)


if __name__ == "__main__":
    # Example usage of TranscriptDesigner
    peptide = "MYPFIRTARMTVCAKKHVHLTRDAAEQLLADIDRRLDQLLPVEGERDVVGAAMREGALAPGKRIRPMLLLLTARDLGCAVSHDGLLDLACAVEMVHAASLILDDMPCMDDAKLRRGRPTIHSHYGEHVAILAAVALLSKAFGVIADADGLTPLAKNRAVSELSNAIGMQGLVQGQFKDLSEGDKPRSAEAILMTNHFKTSTLFCASMQMASIVANASSEARDCLHRFSLDLGQAFQLLDDLTDGMTDTGKDSNQDAGKSTLVNLLGPRAVEERLRQHLQLASEHLSAACQHGHATQHFIQAWFDKKLAAVS"

    designer = TranscriptDesigner()
    designer.initiate()

    ignores = set()
    transcript = designer.run(peptide, ignores)

    cds = ''.join(transcript.codons)
    print(f"CDS length: {len(cds)}")
    print(f"CDS starts with ATG: {cds.startswith('ATG')}")
    print(f"CDS ends with stop:  {cds[-3:] in ('TAA','TGA','TAG')}")
    print(f"RBS: {transcript.rbs.gene_name}")
