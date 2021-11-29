import pathlib 
import math 
import matplotlib.pyplot as plt 

x1, x2, dx, A = -5.12, 5.12, 0.005, 10 


dx *= 1000 
x1 *= 1000 
x2 *= 1000 

x1 += 5120 
x2 += 5120 


x1 = int(x1)     
x2 = int(x2) + 1 
dx = int(dx)     

def f(x, A):     
    return A + x**2 - A*math.cos(2*math.pi*x) 
x = [] 
y = [] 


for i in range(x1,x2, dx):
    v = (i-5120)/1000 
    x.append(v) 
    y.append(f(v,A)) 


    
res = pathlib.Path("results") 
res.mkdir(exist_ok=True) 
file = res / "task1.txt" 

with file.open("w") as f: 
    for a, b in zip(x, y): 
        f.write(f"{a} {b}\n") 

plt.plot(x, y) 
plt.grid() 
plt.savefig("results/task1.png")
