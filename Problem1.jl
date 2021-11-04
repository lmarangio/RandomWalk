using LinearAlgebra

#################  Discretization step  ################### 

k = 10  #<--- distance between two point of the grid, aka "the ant velocity"

x_i = -20
x_f = 20
m = -20
M = 20

n = Int((1 + (x_f-x_i)/k)*(1 + (M-m)/k)) #<--- number of nodes, i.e. the size of the associated transition matrix

#################  Building the matrix  ################### 

A = zeros(n,n)

t = 1 #<--- nodes tracker
origin = 1 #<--- origin tracker
absorbing = 0 #<--- absorbing and transient state counters
transient = 0

for i in 1: Int(1 + (M-m)/k)
    for j in 1 : Int(1 + (x_f-x_i)/k)
        x = x_i + k*(j-1)
        y = M - k*(i-1)

        if x == 0 && y == 0 
            origin = t 
        end
        if x_i<x<x_f && m<y<M  
            A[t, t-1] = 1/4
            A[t, t+1] = 1/4
            A[t,t-Int(1 + (x_f-x_i)/k)] = 1/4
            A[t, t+ Int(1 + (x_f-x_i)/k)] = 1/4
            transient += 1
        else
            for s in 1: n
                A[t,s] = 0
            end
            A[t,t] = 1
            absorbing += 1
        end
        t += 1
    end
end

#################  Canonical form  ################### 
# It is not strictly necesessary to compute the whole canonical form (we just need the transient part); but since we are here...

function swapcol!(x,i,j)
    for k in axes(x, 1)  
        x[k, i], x[k, j] = x[k, j], x[k, i]
    end
end
function swaprow!(x,i,j)
    for k in axes(x, 1) 
        x[i, k], x[j, k] = x[j,k], x[i,k]
    end
end

#First permutation : we put the origin in the first place 
swapcol!(A, 1 ,origin)
swaprow!(A, 1, origin)
#Second permutation : we separate transient state from absorbing state
boo = 1
for i in 1: n
    for j in 1 :n 
        if A[i,j] == 1
            for k in i+1 :n
                for s in 1:n
                    if (A[k,s] == 1/4) && boo == 1
                        swaprow!(A, i, k)
                        swapcol!(A, i, k)
                        boo = 0
                    end
                end
            end
        end
    end
    boo = 1
end

#Third and last permutation (we do not really need this step): rearranging the absorbing matrix into the identity.
B = A[transient+1:n,transient+1:n] 
C = transpose(B)
Z = zeros(n,n)
for i in 1: transient
    Z[i,i] = 1
end
for i in transient +1 : n
    for j in transient + 1 : n
        Z[i,j] = C[i - transient, j - transient]
    end
end
W = Z^-1 * A * Z

#################  Computing the resolvent  ################### 

Q = A[1:transient, 1:transient]
N = (Matrix(1.0I, transient, transient) - Q)^-1 #<---- Resolvent matrix

#Finally, the expected number of steps before being absorbed when starting in transient state i, 
#is the sum of the i-th row of N. We start from the origin, and we permuted the matrix 
#so that the origin is in the first node, hence:

solution = sum(N[1,:])
