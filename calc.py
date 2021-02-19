import numpy as np

INF = 1000 * 1000 * 1000
crible = 1000001 * [INF]

for p in range(2, 1000001):
   if crible[p] == INF:
      for j in range(p, 1000001, p):
         crible[j] = min(crible[j], p)

print("crible ok")

def puis(n, k):
   p = 1
   for i in range(k):
      p *= n
   return p

def phi(n):
   p = 1

   while n != 1:
      d = crible[n]
      k = 0
      
      while n % d == 0:
         n = n // d
         k += 1
      
      p *= puis(d, k) - puis(d, (k - 1))
   
   return p

def sigma(n):
   p = 1
   
   while n != 1:
      k = 0
      d = crible[n]
      
      while n % d == 0:
         n = n // d
         k += 1
      
      p *= (puis(d, (k + 1)) - 1) // (d - 1)
   
   return p

def L(a, b, s):
   return np.sum([puis(phi(i), a) * puis(sigma(i), b) / puis(i, s) for i in range(1, 1000000)])

print(L(0,2,4))
print(L(0,0,2) * (L(0,0,3) ** 2) * L(0,0,4) / L(0,0,6))
