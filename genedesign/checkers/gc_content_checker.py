class GCContentChecker:
    """
    Checks that the GC content of a CDS is within acceptable bounds globally
    (40-65%) and locally in every 50 bp window (30-70%).
    """

    def __init__(self):
        self.global_min = 0.40
        self.global_max = 0.65
        self.local_min = 0.30
        self.local_max = 0.70
        self.window = 50

    def initiate(self):
        pass  # no external data needed

    def run(self, cds: str):
        """
        Returns (passed: bool, gc: float, message: str).
        """
        if not cds:
            return True, 0.0, "empty sequence"

        seq = cds.upper()
        gc = (seq.count('G') + seq.count('C')) / len(seq)

        if not (self.global_min <= gc <= self.global_max):
            msg = f"global GC {gc:.2%} outside [{self.global_min:.0%}, {self.global_max:.0%}]"
            return False, gc, msg

        for i in range(0, len(seq) - self.window + 1, self.window // 2):
            chunk = seq[i:i + self.window]
            local_gc = (chunk.count('G') + chunk.count('C')) / len(chunk)
            if not (self.local_min <= local_gc <= self.local_max):
                msg = f"local GC {local_gc:.2%} at position {i} outside [{self.local_min:.0%}, {self.local_max:.0%}]"
                return False, gc, msg

        return True, gc, "ok"
