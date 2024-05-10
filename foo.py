from gurobipy import *
import numpy as np
import matplotlib.pyplot as plt
from itertools import product

# ST
N = 6

grid_shape = 10, 10
total_vertex_count = grid_shape[0] * grid_shape[1]
grid_size = 10

# x.y 좌표
np.random.seed(1004)
X = list(np.random.randint(low=0, high=grid_size, size=N))
# print(X)
Y = list(np.random.randint(low=0, high=grid_size, size=N))


# src, sink 말고 Grid의 좌표값들

grid = list(product(list(range(grid_size)), repeat=2))
grid = np.array(grid)
# print(grid)
# print(grid.shape)
# print(grid[1, 1])



def draw_ST(X, Y, flow):

    plt.figure(figsize=(12, 5))
    plt.scatter(X, Y, color='blue')

    for i in range(len(X)):
        plt.annotate('$d_{%d}=%d$' % (i, flow[i]), (X[i] - 0.5, Y[i] - 5))

    plt.xlabel("X axis")
    plt.ylabel("Y axis")
    plt.title("nodes")
    plt.show()

# 유량
flow = list(np.random.randint(low=1, high=5, size=N))

# 0 ~ 19 까지 i 값들
ST_index = [i for i in range(N)]
# print(ST_index)
# 0 ~ 99 까지 j 값들
vertex_index = [j for j in range(total_vertex_count)]

# i로부터 j 로 이어질 수 있는 경우의 쌍들
position = [(i, j) for i in ST_index for j in vertex_index]
# print(position)
# 슈퍼 노드의 수: 가정
P = 2

draw_ST(X, Y, flow)
# 거리 행렬
distance_matrix = {(i, j): np.hypot(X[i] - grid[j, 0], Y[i] - grid[j, 1]) for i in ST_index for j in vertex_index}
print(distance_matrix)
print()
print(np.array(distance_matrix).shape)
# P-median
model = Model('P-Median')

# 결정변수

# z ij 이진 변수화 ST의 i 에서 vertex의 j로 이어져 있다면 1 아니면 0
z = model.addVars(position, vtype=GRB.BINARY, name='X')

# y j 이진 변수화 vertex 의 j 번째가 활성화 되어있다면 1 아니면 0
y = model.addVars(vertex_index, vtype=GRB.BINARY, name='Y')

# 목적함수
model.setObjective(quicksum(distance_matrix[i, j] * z[i, j] for i, j in position), GRB.MINIMIZE)

# 제약조건
model.addConstrs(quicksum(z[i, j] for j in vertex_index) == 1 for i in ST_index)
model.addConstr(quicksum(y[j] for j in vertex_index) <= P)
model.addConstrs(z[i, j] - y[j] <= 0 for i in ST_index for j in vertex_index)
model.optimize()


supernodesInVertex = [k for k in vertex_index if y[k].x == 1]
print(supernodesInVertex)


node_to_supernode = [k for k in position if z[k].x == 1]
print(node_to_supernode)

plt.figure(figsize=(12, 5))
plt.scatter(X, Y, color='blue')

for n in supernodesInVertex:
    plt.scatter(grid[n, 0], grid[n, 1], color='green', marker='D', s=150)

for i in range(len(X)):
    plt.annotate('$d_{%d}=%d$' % (i, flow[i]), (X[i] - 0.5, Y[i] - 5))

for n in node_to_supernode:
    i = n[0]
    j = n[1]

    plt.plot([X[i], grid[j, 0]], [Y[i], grid[j, 1]])

plt.xlabel("X axis")
plt.ylabel("Y axis")
# plt.xlim(10)
# plt.ylim(10)
plt.title("Nodes")
plt.grid()
plt.show()
