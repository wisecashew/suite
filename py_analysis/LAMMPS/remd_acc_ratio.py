#!/home/satyend/.conda/envs/phase/bin/python

def parse_remd_output(file_path):
	"""
	Parse the REMD output file and extract replica positions at each step.
	
	:param file_path: Path to the REMD output file
	:return: Dictionary where keys are steps and values are lists of replica positions
	"""
	remd_steps = {}
	
	with open(file_path, 'r') as file:
		for line in file:
			# Identify lines that start with a step number
			parts = line.split()
			if parts and parts[0].isdigit():
				step = int(parts[0])
				replicas = [int(rep) for rep in parts[1:]]
				remd_steps[step] = replicas
				
	return remd_steps


def calculate_acceptance_rate(remd_steps):
	"""
	Calculate the acceptance rate from REMD simulation steps.
	
	:param remd_steps: Dictionary of steps and replica configurations
	:return: Acceptance rate as a float
	"""
	total_attempts = 0
	successful_swaps = 0
	
	# Sort steps to ensure proper order
	steps = sorted(remd_steps.keys())
	
	# Iterate over each consecutive pair of steps
	for i in range(1, len(steps)):
		step_prev = remd_steps[steps[i - 1]]  # Replica configuration at previous step
		step_curr = remd_steps[steps[i]]	  # Replica configuration at current step
		
		# Count how many neighboring pairs have swapped between these steps
		for j in range(len(step_prev) - 1):  # Only consider neighboring pairs (n, n+1)
			total_attempts += 1
			if step_prev[j] != step_curr[j] and step_prev[j+1] != step_curr[j+1]:
				successful_swaps += 1
	
	# Calculate acceptance rate
	print(f"success = {successful_swaps}")
	print(f"total   = {total_attempts}")
	acceptance_rate = successful_swaps / total_attempts if total_attempts > 0 else 0
	return acceptance_rate


if __name__ == "__main__":
	# Example REMD output file path
	file_path = "output"
	
	# Parse REMD output file
	remd_steps = parse_remd_output(file_path)
	
	# Calculate acceptance rate
	acceptance_rate = calculate_acceptance_rate(remd_steps)
	
	# Print the acceptance rate
	print(f"Acceptance Rate: {acceptance_rate:.2}")

