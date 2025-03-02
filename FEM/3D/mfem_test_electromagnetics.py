import mfem.ser as mfem
import numpy as np

# 空気領域のメッシュを作成
mesh = mfem.Mesh(1, 1, 1, "TETRAHEDRON")

# 磁石の位置とサイズを定義
magnet_center = [0.5, 0.5, 0.5]
magnet_size = [0.3, 0.1, 0.05]

# 磁石の領域を定義
magnet_region = mfem.Vector(magnet_size)
magnet_region.Set(1.0, mfem.Vector([1, 1, 1]))

# 有限要素空間を定義
fec = mfem.H1_FECollection(1, mesh.Dimension())   # H1 order=1
fespace = mfem.FiniteElementSpace(mesh, fec)

# 磁界の係数を定義
mu0 = 4 * np.pi * 1e-7  # 真空の透磁率
H_magnet = mfem.VectorConstantCoefficient([0.0, 0.0, 1.0 / mu0])  # 磁石の磁界

# 双線形形式と線形形式を定義
a = mfem.BilinearForm(fespace)
a.AddDomainIntegrator(mfem.DiffusionIntegrator(mfem.ConstantCoefficient(1.0 / mu0)))
a.Assemble()

b = mfem.LinearForm(fespace)
b.AddDomainIntegrator(mfem.DomainLFIntegrator(H_magnet, 2, 0))
b.Assemble()

# 解ベクトルを初期化
x = mfem.GridFunction(fespace)
x.Assign(0.0)

# 線形方程式系を形成 (AX=B)
A = mfem.OperatorPtr()
B = mfem.Vector()
X = mfem.Vector()
a.FormLinearSystem(mfem.intArray(), x, b, A, X, B)
print("Size of linear system: " + str(A.Height()))

# PCGを使って線形方程式を解き、解をxに格納
AA = mfem.OperatorHandle2SparseMatrix(A)
M = mfem.GSSmoother(AA)
mfem.PCG(AA, M, B, X, 1, 200, 1e-12, 0.0)
a.RecoverFEMSolution(X, b, x)

# 頂点と解をnumpy配列として抽出
verts = mesh.GetVertexArray()
sol = x.GetDataArray()

# matplotlibを使って解をプロット
import matplotlib.pyplot as plt
import matplotlib.tri as tri

verts = np.array(verts)
triang = tri.Triangulation(verts[:,0], verts[:,1])

fig, ax = plt.subplots()
ax.set_aspect('equal')
tpc = ax.tripcolor(triang, sol, shading='gouraud')
fig.colorbar(tpc)
plt.show()
