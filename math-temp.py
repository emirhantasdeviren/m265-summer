def MatrixDim(A: list) -> list:
    row = len(A)

    if not A:
        column = 0
    else:
        column = len(A[0])

    return [row, column]


def MatrixMult(A: list, B: list):
    result = []

    dim_a = MatrixDim(A)
    dim_b = MatrixDim(B)

    if dim_a[1] is dim_b[0]:           # Checking if A and B can mutiply
        for i in range(dim_a[0]):
            result.append([])
            for j in range(dim_b[1]):
                c = 0
                for k in range(dim_a[1]):
                    c += A[i][k] * B[k][j]
                result[i].append(c)
        return result
    else:
        print("ERROR: Matrices can not be multiplied")
        return -1


def pivots(A: list) -> list:
    result = []

    for i, row in enumerate(A):
        for j, entry in enumerate(row):
            if entry != 0:
                result.append(j)
                break
        if len(result) < i + 1:
            result.append(-1)
    return result


def M265isRref(A: list) -> bool:
    pivot_points: list = pivots(A)
    result: bool = all(a < b if a != -1 and b != -1 else a >=
                       b for a, b in zip(pivot_points, pivot_points[1:]))

    if result:
        for row, pivot in zip(A, pivot_points):
            for j, entry in enumerate(row):
                if pivot != -1:
                    if j < pivot and entry != 0:
                        return False
                    elif j == pivot and entry != 1:
                        return False
                    elif j > pivot and entry != 0 and j in pivot_points:
                        return False

    return result


def M265BackSubs(A: list, b: list) -> list:
    [row, column] = MatrixDim(A)
    pivot_points: list = pivots(A)
    par_sol: list = [0 for _ in range(column)]
    homg_sol: list = [[0 for _ in range(column)] for i in range(
        column) if i not in pivot_points]
    sol_set: list = []

    if M265isRref(A):
        homg_sol_index: int = 0
        par_sol_index: int = 0
        for j in range(column):
            for i in range(row):
                if j in pivot_points and A[i][j] == 1:
                    par_sol[j] = b[par_sol_index]
                elif pivot_points[i] == -1 and b[i] != 0:
                    return [[], homg_sol]
                elif j not in pivot_points and pivot_points[i] != -1:
                    homg_sol[homg_sol_index][pivot_points[i]] = -A[i][j]
                    homg_sol[homg_sol_index][j] = 1
            if j in pivot_points:
                par_sol_index += 1
            else:
                homg_sol_index += 1
    else:
        print("ERROR: Given matrix is not in reduced row echelon form")

    sol_set.append(par_sol)
    sol_set.append(homg_sol)
    return sol_set


def print_matrix(A: list):
    for i in range(len(A)):
        print(A[i])


def M265Identity(length: int) -> list:
    A = []

    for i in range(length):
        A.append([])
        for j in range(length):
            if i == j:
                A[i].append(1)
            else:
                A[i].append(0)

    return A


def M265GaussianElimination(A: list, debug: bool = False) -> list:
    [row, _column] = MatrixDim(A)
    result = []

    if debug:
        print("Starting")
        print_matrix(A)
        print()


    while not M265isRref(A):
        pivot_points = pivots(A)
        for index, (a, b) in enumerate(zip(pivot_points, pivot_points[1:])):
            if not debug:
                print(
                    f"index = {index}, a = {a}\nindex = {index + 1}, b = {b}")
            if a == -1 and b != -1:
                i = M265Identity(row)
                for j in range(row):
                    i[index][j], i[index + 1][j] = i[index + 1][j], i[index][j]
                result.append(i)
                A = MatrixMult(i, A)
                # pivot_points[index], pivot_points[index +
                                                  # 1] = pivot_points[index + 1], pivot_points[index]
                if debug:
                    print("Found zero row above non-zero row swapping")
                    print_matrix(A)
                    print()
                break
            elif a > b and a != -1 and b != -1:
                i = M265Identity(row)
                for j in range(row):
                    i[index][j], i[index + 1][j] = i[index + 1][j], i[index][j]
                result.append(i)
                A = MatrixMult(i, A)
                # pivot_points[index], pivot_points[index +
                                                  # 1] = pivot_points[index + 1], pivot_points[index]
                if debug:
                    print("Arranging pivot points swapping two adjacent rows")
                    print_matrix(A)
                    print()
                break
            elif a == b and a != -1 and b != -1:
                i = M265Identity(row)
                multiply_by = -A[index + 1][b] / A[index][a]
                for j in range(row):
                    i[index + 1][j] += i[index][j] * multiply_by
                result.append(i)
                A = MatrixMult(i, A)
                if debug:
                    print(
                        "Eliminating pivot point. Multiplying above row with ratio and adding to below row")
                    print_matrix(A)
                break
            else:
                if (a != -1 and A[index][a] != 1) or (b != -1 and A[index + 1][b] != 1):
                    if a != -1 and A[index][a] != 1:
                        i = M265Identity(row)
                        i[index][index] /= A[index][a]
                        result.append(i)
                        A = MatrixMult(i, A)
                        if debug:
                            print("1- Leading variable is not 1. Dividing by ratio")
                            print_matrix(A)
                            print()
                        break
                    if b != -1 and A[index + 1][b] != 1:
                        i = M265Identity(row)
                        i[index + 1][index + 1] /= A[index + 1][b]
                        result.append(i)
                        A = MatrixMult(i, A)
                        if debug:
                            print("1- Leading variable is not one. Dividing by ratio")
                            print_matrix(A)
                            print()
                        break
                else:
                    if debug:
                        print("Cleaning above the leading variable")
                    for i in range(row):
                        if a != -1 and i < index and A[i][a] != 0:
                            identity = M265Identity(row)
                            multiply_by = -A[i][a] / A[index][a]
                            identity[i][index] = identity[index][index] * \
                                multiply_by
                            result.append(identity)
                            A = MatrixMult(identity, A)
                        if b != -1 and i < index + 1 and A[i][b] != 0:
                            identity = M265Identity(row)
                            multiply_by = -A[i][b] / A[index + 1][b]
                            identity[i][index + 1] = identity[index +
                                                       1][index + 1] * multiply_by
                            result.append(identity)
                            A = MatrixMult(identity, A)
                    if debug:
                        print_matrix(A)
                        print()

    return result


def M265GaussianElimination_tests():
    m1 = [[1, 1, 2], [2, 2, 1]]
    m2 = [[0, 0, 0, 0], [1, 1, 2, 0], [2, 2, 4, 1]]
    m3 = [[0, 1, 2, 0], [0, 2, 4, 1]]
    m4 = [[0, 1], [2, 0], [0, 2], [4, 1]]
    m5 = [[0, 1], [2, 0]]
    m6 = [[0, 0], [0, 0], [0, 0]]

    tests = [m6]

    for i in range(len(tests)):
        for elementary_matrix in M265GaussianElimination(tests[i], True):
            tests[i] = MatrixMult(elementary_matrix, tests[i])

    print("result")
    for test in tests:
        print_matrix(test)
        print()


def M265BackSubs_tests():
    print("\n")
    print(" --- TEST 1 ---")
    BS1 = [[1, 1, 0, 3, 4], [0, 0, 1, 1, 8], [0, 0, 0, 0, 0]]
    b1 = [1, 2, 0]
    b1n = [1, 2, 1]
    print(" consistent")
    Sol1 = M265BackSubs(BS1, b1)
    print(Sol1)
    print(" inconsisten")
    Sol1n = M265BackSubs(BS1, b1n)
    print(Sol1n)

    print("\n")
    print(" --- TEST 2 ---")
    print("#var < #equ")
    BS2n = [[0, 1], [0, 0], [0, 0]]
    b2n = [0, 1, 0]
    Sol2n = M265BackSubs(BS2n, b2n)
    print(" inconsistent")
    print(Sol2n)

    print("\n")
    print(" --- TEST 3 ---")
    print("#var < #equ")
    BS2u = [[1, 0], [0, 1], [0, 0]]
    b2u = [3, 4, 0]
    Sol2u = M265BackSubs(BS2u, b2u)
    print(" consistent: unique")
    print(Sol2u)

    print("\n")
    print(" --- TEST 4 ---")
    print("#var < #equ")
    BS2i = [[0, 1], [0, 0], [0, 0]]
    b2i = [3, 0, 0]
    print_matrix(BS2i)
    Sol2i = M265BackSubs(BS2i, b2i)
    print(" consistent: infinite")
    print(Sol2i)

    print("\n")
    print(" --- TEST 5 ---")
    print("#var = #equ")
    BS3n = [[1, 0, 0], [0, 0, 0], [0, 0, 0]]
    b3n = [1, 1, 0]
    print(" inconsistent")
    Sol3n = M265BackSubs(BS3n, b3n)
    print(Sol3n)

    print("\n")
    print(" --- TEST 6 ---")
    print("#var = #equ")
    BS3i = [[0, 1, 0], [0, 0, 1], [0, 0, 0]]
    b3i = [1, 1, 0]
    print(" consistent: infinite ")
    Sol3i = M265BackSubs(BS3i, b3i)
    print(Sol3i)

    print("\n")
    print(" --- TEST 7 ---")
    print("#var = #equ")
    BS3u = [[1, 0, 0], [0, 1, 0], [0, 0, 1]]
    b3u = [2, 1, 3]
    print(" consistent: unique ")
    Sol3u = M265BackSubs(BS3u, b3u)
    print(Sol3u)

    print("\n")
    print(" --- TEST 8 ---")
    BS4 = [[0, 0, 0]]
    b4 = [0]
    b4n = [1]
    Sol4 = M265BackSubs(BS4, b4)
    print("zero matrix (R<C): consistent")
    print(Sol4)
    Sol4n = M265BackSubs(BS4, b4n)
    print("zero matrix (R<C): inconsistent")
    print(Sol4n)

    print("\n")
    print(" --- TEST 9 ---")
    BS5 = [[0], [0]]
    b5 = [0, 0]
    b5n = [1, 0]
    Sol5 = M265BackSubs(BS5, b5)
    print("zero matrix (R>C): consistent")
    print(Sol5)
    Sol5n = M265BackSubs(BS5, b5n)
    print("zero matrix (R>C): inconsistent")
    print(Sol5n)


def mein_test():
    m7 = [[0, 1, 2, 4], [1, 0, 2, 0], [0, 2, 1, 1]]
    print_matrix(m7)
    print()
    for elementary_matrix in M265GaussianElimination(m7):
        print("elemantary matrix:")
        print_matrix(elementary_matrix)
        print()
        print("next matrix")
        m7 = MatrixMult(elementary_matrix, m7)
        print_matrix(m7)
        print()


if __name__ == "__main__":
    M265GaussianElimination_tests()
