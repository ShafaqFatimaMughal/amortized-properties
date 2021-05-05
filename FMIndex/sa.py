from collections import defaultdict
from itertools import zip_longest, islice


def sort_bucket(str, bucket, bucket_level=0, keylen=6):
    d = defaultdict(list)
    for i in bucket:
        key = str[i+bucket_level : i+bucket_level+keylen]
        d[key].append(i)
    result = []
    for k,v in sorted(d.items()):
        if len(v) > 1:
            result += sort_bucket(str, v, bucket_level + keylen, keylen*2)
        elif v:
            result.append(v[0])
    return result

def suffix_array_manber_myers(str):
    """ Returns suffix array without having to create every suffix in the process.
    Much more memory efficient than naive implementation.
    Major drawback is recursion.
    TODO try to replace recursion with iterative approach. """
    return sort_bucket(str, (i for i in range(len(str))))

def to_int_keys_best(l):
    """
    l: iterable of keys
    returns: a list with integer keys
    """
    seen = set()
    ls = []
    for e in l:
        if not e in seen:
            ls.append(e)
            seen.add(e)
    ls.sort()
    index = {v: i for i, v in enumerate(ls)}
    return [index[v] for v in l]

def inverse_array(l):
    n = len(l)
    ans = [0] * n
    for i in range(n):
        ans[l[i]] = i
    return ans

def suffix_array_best(s):
    """
    suffix array of s
    O(n * log(n)^2)
    """
    n = len(s)
    k = 1
    line = to_int_keys_best(s)
    while max(line) < n - 1:
        line = to_int_keys_best(
            [a * (n + 1) + b + 1
             for (a, b) in
             zip_longest(line, islice(line, k, None),
                         fillvalue=-1)])
        k <<= 1
    return inverse_array(line)