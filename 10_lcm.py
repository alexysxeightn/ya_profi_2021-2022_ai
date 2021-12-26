def divisors(n):
    ans = set()
    for i in range(1, int(n ** 0.5) + 1):
        if not n % i:
            ans.add(i)
            ans.add(n // i)
    return ans


def get_ans(n, m):
    set_n = divisors(n)
    set_m = divisors(m)
    ans = set()
    for i in set_n:
        for j in set_m:
            ans.add(i + j)
    return len(ans)


def main():
    t = int(input())
    for _ in range(t):
        a, b = map(int, input().split())
        print(get_ans(a, b))


if __name__ == '__main__':
    main()
