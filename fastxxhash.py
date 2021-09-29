# Derived from this implementation of XXHASH64: https://github.com/Cyan4973/xxHash/blob/dev/xxhash.h

from xxhash import xxh64

XXH_PRIME64_1 = 0x9E3779B185EBCA87
XXH_PRIME64_2 = 0xC2B2AE3D27D4EB4F
XXH_PRIME64_3 = 0x165667B19E3779F9
XXH_PRIME64_4 = 0x85EBCA77C2B2AE63
XXH_PRIME64_5 = 0x27D4EB2F165667C5

# the multiplicative inverses of the above, mod 1<<64
# calculated using inverse_mod in sagemath
INV_1 = 0x887493432badb37
INV_2 = 0xba79078168d4baf
INV_3 = 0xe9e9f4c41d6df849
INV_4 = 0xd872e78f6fe1434b
INV_5 = 0xc592c09fdfba7f0d

MASK64 = (1<<64)-1

def XXH_rotl64(x, n):
	return ((x << n) | (x >> (64 - n))) & MASK64

def inverse_round(acc, inp):
	inp = (inp * INV_1) & MASK64
	inp = XXH_rotl64(inp, 64-31)
	inp = (inp - acc) & MASK64
	inp = (inp * INV_2) & MASK64
	return inp & MASK64

def XXH64_round(acc, inp):
	acc += (inp * XXH_PRIME64_2)
	acc &= MASK64
	acc  = XXH_rotl64(acc, 31)
	acc *= XXH_PRIME64_1
	acc &= MASK64
	return acc

def XXH64_mergeRound(acc, val):
	val  = XXH64_round(0, val)
	acc ^= val
	acc  = acc * XXH_PRIME64_1 + XXH_PRIME64_4
	acc &= MASK64
	return acc

def inverse_avalanche(h64):
	h64 ^= h64 >> 32
	h64 *= INV_3
	h64 &= MASK64
	h64 ^= (h64 >> 29) ^ (h64 >> 58)
	h64 *= INV_2
	h64 &= MASK64
	h64 ^= h64 >> 33
	return h64

def inverse_finalize32(prefinal, postfinal):
	h64 = postfinal
	h64 = (h64 - XXH_PRIME64_3) & MASK64
	h64 = (h64 * INV_2) & MASK64
	h64 = XXH_rotl64(h64, 64-23)
	h64 = ((h64 ^ prefinal) * INV_1) & MASK64
	return h64

def inverse_finalize64(prefinal, postfinal):
	h64 = postfinal
	h64 = (h64 - XXH_PRIME64_4) & MASK64
	h64 = (h64 * INV_1) & MASK64
	h64 = XXH_rotl64(h64, 64-27)
	h64 = inverse_round(0, h64 ^ prefinal)
	return h64

# TODO: handle longer suffixes
def inverse_suffix(h64, buf):
	for x in buf[::-1]:
		h64  = (h64 * INV_1) & MASK64
		h64  = XXH_rotl64(h64, 64-11)
		h64  ^= (x * XXH_PRIME64_5) & MASK64
	return h64

def premine(data, seed, total_len):
	if not data:
		return (seed + XXH_PRIME64_5 + total_len) & MASK64
	v1 = (seed + XXH_PRIME64_1 + XXH_PRIME64_2) & MASK64
	v2 = (seed + XXH_PRIME64_2) & MASK64
	v3 = seed
	v4 = (seed - XXH_PRIME64_1) & MASK64
	
	for i in range(0, len(data)-31, 32):
		v1 = XXH64_round(v1, int.from_bytes(data[i+0 :][:8], "little"))
		v2 = XXH64_round(v2, int.from_bytes(data[i+8 :][:8], "little"))
		v3 = XXH64_round(v3, int.from_bytes(data[i+16:][:8], "little"))
		v4 = XXH64_round(v4, int.from_bytes(data[i+24:][:8], "little"))
	
	h64 = XXH_rotl64(v1, 1) + XXH_rotl64(v2, 7) + XXH_rotl64(v3, 12) + XXH_rotl64(v4, 18)
	h64 &= MASK64
	
	h64 = XXH64_mergeRound(h64, v1)
	h64 = XXH64_mergeRound(h64, v2)
	h64 = XXH64_mergeRound(h64, v3)
	h64 = XXH64_mergeRound(h64, v4)
	
	h64 += total_len
	h64 &= MASK64
	
	buf = data[i+32:]
	
	while len(buf) >= 8:
		x, buf = buf[:8], buf[8:]
		#print(buf)
		x = int.from_bytes(x, "little")
		h64 ^= XXH64_round(0, x)
		h64  = XXH_rotl64(h64,27) * XXH_PRIME64_1 + XXH_PRIME64_4
		h64 &= MASK64
	
	assert(len(buf) == 0)
	
	return h64


def fastmine_inner(prefix, target, suffix=b"", numeric_only=True, seed=2811, brutelen=8):
	postfinal = inverse_avalanche(target)
	postfinal = inverse_suffix(postfinal, suffix)
	prefinal = premine(prefix, seed, len(prefix)+brutelen+len(suffix))
	if brutelen == 8:
		nonce = inverse_finalize64(prefinal, postfinal)
		nonce = nonce.to_bytes(8, "little")
	elif brutelen == 4:
		nonce = inverse_finalize32(prefinal, postfinal)
		if nonce >= 0x100000000:
			return None
		nonce = nonce.to_bytes(4, "little")
	else:
		assert(False)
	
	# DUCO PoW results must have a numeric suffix
	if numeric_only:
		good = True
		for c in nonce:
			if c < 0x30 or c > 0x39:
				good = False
				break
		if not good:
			return None
	lol = prefix + nonce + suffix

	if numeric_only:
		return int((nonce + suffix).decode())
	else:
		return nonce + suffix


def fastmine(prefix, target):
	nonce = fastmine_inner(prefix, target)
	if nonce is not None:
		return nonce
	
	for i in range(10):
		nonce = fastmine_inner(prefix, target, str(i).encode())
		if nonce is not None:
			return nonce
	
	for j in [3, 2, 1]:
		for i in range(10**j):
			nonce = fastmine_inner(prefix, target, str(i).rjust(j, "0").encode(), brutelen=4)
			if nonce is not None:
				return nonce
	
	nonce = fastmine_inner(prefix, target, brutelen=4)
	if nonce is not None:
		return nonce
	
	for i in range(1000):
		if xxh64(prefix + str(i).encode(), seed=2811).intdigest() == target:
			return i
	
	return None



if __name__ == "__main__":
	print("RUNNING TESTS")
	print()

	target_hash = xxh64(b"PREFIXES"*5 + b"123456789", seed=2811).intdigest()
	print("hash inversion:")
	print(fastmine(b"PREFIXES"*5, target_hash))
	print()

	target_hash = xxh64(b"PREFIXES"*5 + b"12345678", seed=2811).intdigest()
	print("hash inversion:")
	print(fastmine(b"PREFIXES"*5, target_hash))
	print()

	target_hash = xxh64(b"PREFIXES"*5 + b"1234567", seed=2811).intdigest()
	print("hash inversion:")
	print(fastmine(b"PREFIXES"*5, target_hash))
	print()

	target_hash = xxh64(b"PREFIXES"*5 + b"123456", seed=2811).intdigest()
	print("hash inversion:")
	print(fastmine(b"PREFIXES"*5, target_hash))
	print()

	target_hash = xxh64(b"PREFIXES"*5 + b"12345", seed=2811).intdigest()
	print("hash inversion:")
	print(fastmine(b"PREFIXES"*5, target_hash))
	print()

	target_hash = xxh64(b"PREFIXES"*5 + b"1234", seed=2811).intdigest()
	print("hash inversion:")
	print(fastmine(b"PREFIXES"*5, target_hash))
	print()

	target_hash = xxh64(b"PREFIXES"*5 + b"123", seed=2811).intdigest()
	print("hash inversion:")
	print(fastmine(b"PREFIXES"*5, target_hash))
	print()

	target_hash = xxh64(b"PREFIXES"*5 + b"12", seed=2811).intdigest()
	print("hash inversion:")
	print(fastmine(b"PREFIXES"*5, target_hash))
	print()

	target_hash = xxh64(b"PREFIXES"*5 + b"1", seed=2811).intdigest()
	print("hash inversion:")
	print(fastmine(b"PREFIXES"*5, target_hash))
	print()
