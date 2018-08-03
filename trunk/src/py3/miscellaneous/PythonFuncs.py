def logical_and_n(*arg):
	from numpy import logical_and

	if len(arg):
		result=arg[0]
		for item in arg[1:]:
			result=logical_and(result,item)
		return result

	else:
		return None

def logical_or_n(*arg):
	from numpy import logical_or

	if len(arg):
		result=arg[0]
		for item in arg[1:]:
			result=logical_or(result,item)
		return result

	else:
		return None

