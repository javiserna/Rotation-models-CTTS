import matplotlib.pyplot as plt

model1=[0.764,0.838,0.744,0.920]
model2=[0.675,0.478,0.350,0.367]
model3=[0.771,0.434,0.362,0.367]
time=[1,2,3,4]


plt.plot(time, model1, label=r'Model $B_{\ast}$ decreases with age')
plt.plot(time, model2, label=r'Model $B_{\ast}$ increases with age')
plt.plot(time, model3, label=r'Model $B_{\ast}$ constant with age')

plt.ylabel('Similarity')
plt.xlabel('bin')
plt.ylim(0, 1)
plt.legend()
#plt.show()
plt.savefig("Similarity_Bfield.png", dpi=300)
