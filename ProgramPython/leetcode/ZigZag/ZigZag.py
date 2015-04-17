class Solution:
    # @return a string
    def convert(self, s, nRows):
        if nRows == 1:
            return s
        D = 2*nRows -2
        L = len(s)
        R = ''
        for i in range(0, nRows):
            MD = 2*nRows - 2 - 2*i
            t = i
            while t < L:
                R += s[t]
                if i != 0 and i!= nRows-1 and t+MD < L:
                    R += s[t+MD]
                t += D
        return R

sp = Solution()
p = sp.convert('PAYPALISHIRING', 3)
print p
