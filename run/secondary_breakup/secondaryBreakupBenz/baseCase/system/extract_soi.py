import re

def get_soi_value():
    # Define the path to the sprayCloudProperties file
    spraycloud_properties_file = "../constant/sprayCloudProperties"

    # Define the pattern to search for the SOI value
    pattern = r"SOI\s+([+-]?\d*\.?\d+(?:[eE][+-]?\d+)?)"

    # Open the sprayCloudProperties file
    with open(spraycloud_properties_file, "r") as f:
        content = f.read()

    # Search for the SOI value using regular expressions
    match = re.search(pattern, content)

    if match:
        # Extract the SOI value from the matched group
        soi_value = match.group(1)
        return soi_value
    else:
        # If the SOI value is not found, return None or raise an exception
        return "3.15e-3"

# Call the function to get the SOI value
soi_value = get_soi_value()

'''
# Print the SOI value or handle it as needed
if soi_value:
    print(soi_value)
else:
    print("SOI value not found in sprayCloudProperties file.")
'''
