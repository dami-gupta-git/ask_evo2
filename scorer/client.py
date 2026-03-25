import requests

ENDPOINT_URL = "https://dami-gupta-git--askevo2-scorer-score-variant.modal.run"


def score_variant(ref_sequence: str, alt_sequence: str) -> dict:
    response = requests.post(
        ENDPOINT_URL,
        json={"ref_sequence": ref_sequence, "alt_sequence": alt_sequence},
        timeout=300,
    )
    response.raise_for_status()
    return response.json()

__all__ = ["score_variant"]
