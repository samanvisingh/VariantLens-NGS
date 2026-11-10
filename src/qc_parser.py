"""
VariantLens-NGS: FASTQ Quality Control & Reader Module
Author: Samanvi Singh
"""

import os
from typing import Dict, List, Tuple


class FASTQReader:
    """Parses and computes quality statistics for raw NGS FASTQ reads."""

    def __init__(self, filepath: str):
        self.filepath = filepath
        if not os.path.exists(filepath):
            raise FileNotFoundError(f"FASTQ file not found at: {filepath}")

    @staticmethod
    def phred_score(char: str) -> int:
        """Converts an ASCII character to a Phred+33 quality score.

        Formula: Q = ord(char) - 33
        """
        return ord(char) - 33

    def parse_reads(self) -> List[Dict[str, str]]:
        """Parses FASTQ file into a list of structured read dictionaries.

        FASTQ format structure per read:
        Line 1: @Header
        Line 2: Sequence String
        Line 3: +
        Line 4: Quality ASCII String
        """
        reads = []
        with open(self.filepath, "r") as handle:
            while True:
                header = handle.readline().strip()
                if not header:
                    break  # End of file reached
                sequence = handle.readline().strip()
                plus_line = handle.readline().strip()
                quality = handle.readline().strip()

                # Basic validation for standard FASTQ integrity
                if not header.startswith("@"):
                    raise ValueError(f"Invalid header line: {header}")
                if len(sequence) != len(quality):
                    raise ValueError(
                        f"Sequence length ({len(sequence)}) does not match quality length ({len(quality)}) for read {header}"
                    )

                reads.append(
                    {
                        "header": header,
                        "sequence": sequence,
                        "quality_str": quality,
                        "phred_scores": [self.phred_score(c) for c in quality],
                    }
                )
        return reads

    def compute_summary_qc(self, min_q_threshold: int = 20) -> Dict[str, float]:
        """Calculates core Quality Control metrics across all reads in the sample.

        - Q20/Q30 percentages: Ratio of bases meeting minimum quality benchmarks.
        - Average Read Quality: Overall mean Phred score across reads.
        """
        reads = self.parse_reads()
        if not reads:
            return {"total_reads": 0, "q20_rate": 0.0, "q30_rate": 0.0, "mean_phred": 0.0}

        total_bases = 0
        q20_bases = 0
        q30_bases = 0
        phred_sum = 0

        for read in reads:
            scores = read["phred_scores"]
            total_bases += len(scores)
            phred_sum += sum(scores)
            q20_bases += sum(1 for q in scores if q >= 20)
            q30_bases += sum(1 for q in scores if q >= 30)

        return {
            "total_reads": len(reads),
            "total_bases": total_bases,
            "q20_rate": round((q20_bases / total_bases) * 100, 2),
            "q30_rate": round((q30_bases / total_bases) * 100, 2),
            "mean_phred": round(phred_sum / total_bases, 2),
        }
