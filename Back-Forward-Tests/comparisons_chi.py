import matplotlib.pyplot as plt

model1=[0.728,0.871,0.915,0.890]
model2=[0.701,0.471,0.350,0.334]
model3=[0.795,0.434,0.362,0.367]
time=[1,2,3,4]


plt.plot(time, model1, label=r'Model $\chi_{\ast}$ decreases with age')
plt.plot(time, model2, label=r'Model $\chi_{\ast}$ increases with age')
plt.plot(time, model3, label=r'Model $\chi_{\ast}$ constant with age')

plt.ylabel('Similarity')
plt.xlabel('bin')
plt.ylim(0, 1)
plt.legend()
#plt.show()
plt.savefig("Similarity_chi.png", dpi=300)
