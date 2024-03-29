import random
import numpy as np
import sys


if len(sys.argv) != 3:
    print("Usage: python3 FindFeasSol.py m n")
    sys.exit(1)
m = int(sys.argv[1])  # Number of rows
n = int(sys.argv[2])  # Number of columns
# random.seed(42)
A = np.array( [[random.gauss(0, 1) for _ in range(n)] for _ in range(m)] )
c0 = np.array( [random.gauss(0, 1) for _ in range(n)] )
magnitude = np.linalg.norm(c0)
c = c0 / magnitude
# print(c)
x0 = ( 1.0/np.sqrt(2*np.log((m*1.0)/(n*1.0))) ) * c
x0_norm = np.linalg.norm(x0)
print("x0_norm =", x0_norm)


def findViolations(x,epsilon):
    b = np.dot(A, x)
    violated_idx=[]
    for i in range(m):
        if b[i] > 1-epsilon:
            violated_idx.append(i)
    return b , violated_idx

def generateNewVector(b, epsilon):
    b1 = [] 
    M0 =  [] 
    for i in violated_idx:
        b1.append( 1.0 - epsilon - b[i])
        M0.append( A[i]- (np.dot(A[i], c)) * c )
    M=np.transpose( np.array(M0)) 
    b1 = np.array(b1)
    u = np.linalg.solve( np.dot(np.transpose(M), M) , b1)
    x1 = np.dot(M,u)
    return x1

max_iter = 30
i = 0
epsilon = 0.1
sum_x0_till_xi = x0
b , violated_idx = findViolations(x0,epsilon)
xi_norm = np.linalg.norm(x0)
num_vc=[] # num of violated constraints
while len(violated_idx) > 0 and i < max_iter: 
    print(f"sum_x0_till_x{i} is not feasible.")
    num_vc.append( len(violated_idx))
    xi = generateNewVector(b,epsilon)
    xi_norm = np.linalg.norm(xi)
    sum_x0_till_xi = [s+h for s,h in zip(sum_x0_till_xi , xi)]
    epsilon = epsilon * 0.1
    b , violated_idx = findViolations(sum_x0_till_xi, epsilon)
    num_vc.append( len(violated_idx))
    i += 1

if len(violated_idx) == 0:
    obj = np.dot(sum_x0_till_xi , c) 
    print(f"i = {i} obj = {obj:.6f}")
    file_path = "results_table.txt"
    with open(file_path , 'a') as file: # append mode
        file.write(f"{m} & {n} & {i} & {x0_norm:.6f} & {num_vc[0]} & {num_vc[1]} \\\ \n")
print("end of while loop")




