import numpy as np
import matplotlib.pyplot as plt

# Problem 7.9______________________________________________________________________

# Part A


rows = []

with open('blur.txt','r') as file:
    for line in file:
        if line.strip():
            row = [float(value) for value in line.split()]
            rows.append(row) 

data = np.array(rows)


print(data.shape)


plt.imshow(data, cmap='gray')


plt.xlabel('x (pixels)')
plt.ylabel('y (pixels)')
plt.title('Blurred Photo')
plt.colorbar(label='Brightness')


plt.show()

#___________________________________________________________________________________________________


# Part B

# f(x,y) = exp(- (x^2+y^2)/(2(sigma)^2)) ; sigma = 25 

sigma = 25
rows_num = data.shape[0]
cols_num = data.shape[1]

psf = np.zeros((rows_num, cols_num))

for i in range(rows_num):
    for j in range(cols_num):

        if i < rows_num/2:
            y = i
        else:
            y = i - rows_num

        if j < cols_num/2:
            x = j
        else:
            x = j - cols_num

        psf[i, j] = np.exp(-(x**2 + y**2) / (2*sigma**2))



plt.imshow(psf, cmap='gray')


plt.xlabel('x (pixels)')
plt.ylabel('y (pixels)')
plt.title('Point Spread Function')
plt.colorbar(label='Brightness')


plt.show()


#____________________________________________________________________________________________________


# Part C

# FT(blurred) = FT(true photo) x FT(PSF)
# FT(true photo) = FT(blurred)/FT(PSF)
# then do inverse transform 


epsilon = 1e-3

blur_fft = np.fft.rfft2(data)
psf_fft = np.fft.rfft2(psf)

unblur_fft = np.zeros(blur_fft.shape, dtype=complex)

for i in range(blur_fft.shape[0]):
    for j in range(blur_fft.shape[1]):

        if abs(psf_fft[i, j]) > epsilon:
            unblur_fft[i, j] = blur_fft[i, j] / psf_fft[i, j]
        else:
            unblur_fft[i, j] = blur_fft[i, j]

unblurred = np.fft.irfft2(unblur_fft, s=data.shape)


plt.imshow(unblurred, cmap='gray')


plt.xlabel('x (pixels)')
plt.ylabel('y (pixels)')
plt.title('Unblurred Photo')
plt.colorbar(label='Brightness')


plt.show()






















