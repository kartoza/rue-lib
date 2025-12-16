# src/rue_lib/financial.py

import json
from pathlib import Path
from typing import Any


class FinancialModel:
    """Financial model class."""

    def data_in_json(self) -> dict[str, Any]:
        """Return data in json."""
        return self.__dict__

    def save(self, output_dir: str):
        """Save financial attributes to the financial.json."""
        output_dir = Path(output_dir)
        output_dir.mkdir(parents=True, exist_ok=True)
        output_path = output_dir / "financial.json"
        output_path.write_text(json.dumps(self.data_in_json(), indent=2))
