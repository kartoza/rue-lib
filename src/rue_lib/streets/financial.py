# src/rue_lib/streets/financial.py

from rue_lib.financial import FinancialModel
from rue_lib.streets import StreetConfig


class FinancialStreet(FinancialModel):
    """Financial attributes of a street."""

    def __init__(self, config=StreetConfig):
        """Initialize a FinancialStreet object."""
        self.save(output_dir=config.output_dir)
