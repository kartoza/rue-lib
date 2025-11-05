"""Basic usage examples for rue-lib."""

from rue_lib.core import display_info, greet

# Simple greeting
print(greet("Python Developer"))
print(greet())  # Default greeting

# Display library information with rich formatting
print("\nLibrary Information:")
display_info()
