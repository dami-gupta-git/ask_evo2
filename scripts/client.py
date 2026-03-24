import requests


ENDPOINT_URL = "https://dami-gupta-git--askevo2-scorer-score-variant.modal.run"


def score_variant(ref_sequence: str, alt_sequence: str) -> dict:
    response = requests.post(
        ENDPOINT_URL,
        json={
            "ref_sequence": ref_sequence,
            "alt_sequence": alt_sequence,
        },
        timeout=300,
    )w
    response.raise_for_status()
    return response.json()


if __name__ == "__main__":
    BRCA1_REF = "ATGGATTTATCTGCTCTTCGCGTTGAAGAAGTACAAAATGTCATTAATGCTATGCAGAAA"
    BRCA1_ALT = "ATGGATTTATCTGCTCTTCGCGTTGAAGAAGTACAAAATGTCATTAATGCTATGCAGAAG"

    result = score_variant(BRCA1_REF, BRCA1_ALT)
    print(f"ref_ll:         {result['ref_ll']:.4f}")
    print(f"alt_ll:         {result['alt_ll']:.4f}")
    print(f"delta:          {result['delta']:.4f}")
    print(f"interpretation: {result['interpretation']}")