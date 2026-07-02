"""Phase 8.2: Literature reference data loader.

Per-target known-planet tables and pipeline config live in
data/known_planets.json so that adding a survey target means editing
data, not code. See the _schema entry in the JSON for field meaning.
"""

import json
import logging
from pathlib import Path

logger = logging.getLogger(__name__)

DATA_PATH = (Path(__file__).resolve().parent.parent
             / "data" / "known_planets.json")

_cache = None


def _load():
    global _cache
    if _cache is None:
        with open(DATA_PATH, encoding="utf-8") as f:
            data = json.load(f)
        data.pop("_schema", None)
        for hd, entry in data.items():
            planets = entry.get("planets", [])
            for p in planets:
                if "period_days" not in p:
                    raise ValueError(
                        f"known_planets.json: planet without period_days "
                        f"under '{hd}'"
                    )
        _cache = data
        logger.info("Loaded reference data for %d targets from %s",
                    len(data), DATA_PATH.name)
    return _cache


def load_known_planets(hd_id):
    """Return the literature planet list for an HD identifier, or None.

    Each planet dict has at least period_days; k_ms, eccentricity,
    confirmed and reference drive fit initialization and
    disputed-signal rejection when present.
    """
    entry = _load().get(hd_id)
    if entry is None:
        return None
    return entry.get("planets") or None


def load_target_config(hd_id):
    """Return per-target pipeline config (possibly empty dict).

    Recognized keys: gp_hyperparams, exclude_instruments,
    activity_cycle_days, rotation_days, notes.
    """
    entry = _load().get(hd_id)
    if entry is None:
        return {}
    return entry.get("config", {})


def known_target_ids():
    """All HD identifiers with reference data."""
    return sorted(_load().keys())
