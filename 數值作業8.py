# 題目資料
x = [4.0, 4.2, 4.5, 4.7, 5.1, 5.5, 5.9, 6.3]
y = [102.6, 113.2, 130.1, 142.1, 167.5, 195.1, 224.9, 256.8]
n = len(x)

### 幫助函數 ###
def solve_2x2(A, B):
    a, b = A[0]
    c, d = A[1]
    det = a * d - b * c
    a1 = (B[0] * d - b * B[1]) / det
    a2 = (a * B[1] - B[0] * c) / det
    return a1, a2

def solve_3x3(A, B):
    def det3(m):
        return (m[0][0]*(m[1][1]*m[2][2] - m[1][2]*m[2][1]) -
                m[0][1]*(m[1][0]*m[2][2] - m[1][2]*m[2][0]) +
                m[0][2]*(m[1][0]*m[2][1] - m[1][1]*m[2][0]))
    D = det3(A)
    def col_replace(mat, col, new_col):
        return [
            [new_col[0] if i == col else mat[0][i] for i in range(3)],
            [new_col[1] if i == col else mat[1][i] for i in range(3)],
            [new_col[2] if i == col else mat[2][i] for i in range(3)]
        ]
    D0 = det3(col_replace(A, 0, B))
    D1 = det3(col_replace(A, 1, B))
    D2 = det3(col_replace(A, 2, B))
    return D0 / D, D1 / D, D2 / D

# === 1a: Least squares polynomial ===
Sx = Sx2 = Sx3 = Sx4 = Sy = Sxy = Sx2y = 0
for i in range(n):
    xi = x[i]
    yi = y[i]
    xi2 = xi * xi
    xi3 = xi2 * xi
    xi4 = xi3 * xi
    Sx += xi
    Sx2 += xi2
    Sx3 += xi3
    Sx4 += xi4
    Sy += yi
    Sxy += xi * yi
    Sx2y += xi2 * yi
A = [[n, Sx, Sx2], [Sx, Sx2, Sx3], [Sx2, Sx3, Sx4]]
B = [Sy, Sxy, Sx2y]
a0, a1, a2 = solve_3x3(A, B)
print("1a. p(x) = {:.6f} + {:.6f}x + {:.6f}x^2".format(a0, a1, a2))
error = sum((y[i] - (a0 + a1 * x[i] + a2 * x[i]**2))**2 for i in range(n))
print("1a Error = {:.6f}".format(error))

# === 1b: f(x) = b * e^(a*x) ===
Sx = Sy = Sxy = Sx2 = 0
lny = []
for i in range(n):
    lny.append(0)
    val = y[i]
    while val > 2.718281828459045:
        val /= 2.718281828459045
        lny[i] += 1
    while val < 1:
        val *= 2.718281828459045
        lny[i] -= 1
    term = val - 1
    lny[i] += term - term**2/2 + term**3/3 - term**4/4
    Sx += x[i]
    Sy += lny[i]
    Sxy += x[i] * lny[i]
    Sx2 += x[i]**2
a, ln_b = solve_2x2([[Sx2, Sx], [Sx, n]], [Sxy, Sy])
b = 2.718281828459045 ** ln_b
print("1b. f(x) = {:.6f} * e^({:.6f}x)".format(b, a))
error = sum((y[i] - b * (2.718281828459045 ** (a * x[i])))**2 for i in range(n))
print("1b Error = {:.6f}".format(error))

# === 1c: f(x) = b * x^a ===
Sx = Sy = Sxy = Sx2 = 0
lnx, lny = [], []
for i in range(n):
    # ln(x)
    lx = 0
    val = x[i]
    while val > 2.718281828459045:
        val /= 2.718281828459045
        lx += 1
    while val < 1:
        val *= 2.718281828459045
        lx -= 1
    tx = val - 1
    lx += tx - tx**2/2 + tx**3/3 - tx**4/4
    # ln(y)
    ly = 0
    val = y[i]
    while val > 2.718281828459045:
        val /= 2.718281828459045
        ly += 1
    while val < 1:
        val *= 2.718281828459045
        ly -= 1
    ty = val - 1
    ly += ty - ty**2/2 + ty**3/3 - ty**4/4
    lnx.append(lx)
    lny.append(ly)
    Sx += lx
    Sy += ly
    Sx2 += lx * lx
    Sxy += lx * ly
a, ln_b = solve_2x2([[Sx2, Sx], [Sx, n]], [Sxy, Sy])
b = 2.718281828459045 ** ln_b
print("1c. f(x) = {:.6f} * x^{:.6f}".format(b, a))
error = sum((y[i] - b * (x[i] ** a))**2 for i in range(n))
print("1c Error = {:.6f}".format(error))

# === 2: Degree-2 Polynomial Least Squares on [-1,1] ===
def f2(x):
    t2 = x*x
    t4 = t2*t2
    cos = 1 - t2/2 + t4/24
    sin2x = 2*x - 8*x**3/3
    return 0.5*cos + 0.25*sin2x
m = 100
h = 2.0 / m
sx0 = sx1 = sx2 = sx3 = sx4 = sf = sfx = sfx2 = 0
for i in range(m + 1):
    xi = -1 + h * i
    fi = f2(xi)
    sx0 += 1
    sx1 += xi
    sx2 += xi * xi
    sx3 += xi**3
    sx4 += xi**4
    sf += fi
    sfx += fi * xi
    sfx2 += fi * xi * xi
a0, a1, a2 = solve_3x3([[sx0, sx1, sx2], [sx1, sx2, sx3], [sx2, sx3, sx4]], [sf, sfx, sfx2])
print("2. Least squares polynomial on [-1,1]: {:.6f} + {:.6f}x + {:.6f}x^2".format(a0, a1, a2))

# === 3: Trigonometric Least Squares for f(x)=x^2*sin(x) on [0,1] ===
def f3(x): return x*x*(x - x**3/6)
m = 16
pi = 3.141592653589793
s = [0]*9
for j in range(m):
    xj = j / m
    fj = f3(xj)
    s[0] += fj
    for k in range(1, 5):
        angle = 2 * pi * k * xj
        sin = angle - (angle**3)/6 + (angle**5)/120
        cos = 1 - (angle**2)/2 + (angle**4)/24
        s[2*k-1] += fj * cos
        s[2*k] += fj * sin
a0 = s[0] / m
a = [s[i] * 2 / m for i in range(1, 9)]
def S4(x):
    result = a0
    for k in range(1, 5):
        angle = 2 * pi * k * x
        sin = angle - (angle**3)/6 + (angle**5)/120
        cos = 1 - (angle**2)/2 + (angle**4)/24
        result += a[2*k-2] * cos + a[2*k-1] * sin
    return result
I1 = I2 = 0
steps = 100
for i in range(steps + 1):
    xi = i / steps
    w = 1
    if i == 0 or i == steps:
        w = 0.5
    I1 += w * S4(xi)
    I2 += w * f3(xi)
I1 /= steps
I2 /= steps
print("3b. ∫S4(x)dx on [0,1]: {:.6f}".format(I1))
print("3c. ∫x²sin(x)dx on [0,1]: {:.6f}".format(I2))
print("3d. E(S4) = {:.6f}".format((I2 - I1)**2))
