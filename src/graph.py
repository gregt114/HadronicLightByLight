import numpy as np
import matplotlib.pyplot as plt


lmdv = np.loadtxt("dataLMDV.txt", delimiter=",")
q4   = np.loadtxt("dataQ4.txt", delimiter=",")
q6   = np.loadtxt("dataQ6.txt", delimiter=",")


plt.figure(figsize=[10,8])


plt.plot(lmdv[:,0], lmdv[:,1], "k-", label="LMDV")
plt.plot(q4[:,0], q4[:,1], "r-", label="Q4")
plt.plot(q6[:,0], q6[:,1], "b-", label="Q6")

plt.title("Form Factor Values")
plt.xlabel("Q1 = Q2")
plt.ylabel("Form Factor Value")

plt.grid()
plt.legend()
plt.show()
