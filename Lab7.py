import re

# Original string
text = '$17£29€8¥55'

# Regex with two groups:
# Group 1 = currency symbol
# Group 2 = numeric value
pattern = r'([$£€¥])(\d+)'

# Exchange rates to British pounds
# Dictionary: exchange_rates stores currency symbols as keys and exchange rates as values.
exchange_rates = {
    '$': 0.82,
    '£': 1.00,
    '€': 0.87,
    '¥': 0.11
}

# Variable to store the total value in pounds
total_pounds = 0

# Find all matches
# Loop: for match in re.finditer(...) processes each match.
for match in re.finditer(pattern, text):

    # Extract the two groups
    # Capturing groups: group(1) and group(2) extract different parts of the match.
    # Integer conversion: int() converts the extracted number from a string to an integer.
    currency = match.group(1)
    value = int(match.group(2))

    # Convert to British pounds
    pounds = value * exchange_rates[currency]

    # Add to total
    # Accumulation: total_pounds += pounds adds each converted value to the running total.
    total_pounds += pounds

    # Print each conversion
    # Formatting: :.2f displays a number to two decimal places.
    print(f"{currency}{value} = £{pounds:.2f}")

# Print the total
print(f"Total value in pounds = £{total_pounds:.2f}")
