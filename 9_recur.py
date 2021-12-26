INF = 10 ** 9 + 7

def tr(A):
  return A[0][0] + A[1][1]

def det(A):
  return A[0][0] * A[1][1] - A[0][1] * A[1][0] 

def determinant(a, b, c):
  return b ** 2 - 4 * a * c

def find_type_recur(A):
  tr_A, det_A = tr(A), det(A)
  D = determinant(1, -tr_A, det_A)
  if D > 0:
    return 0 # lambda_1 != lambda_2, lambda Real
  elif 0 == D:
    return 1 # lambda_1 == lambda_2, lambda Real
  return 2 # lambda Complex

def eigenvectors_and_eigenvalues_first_type(A):
  tr_A, det_A = tr(A), det(A)
  D = determinant(1, -tr_A, det_A)

  lambda_1, lambda_2 = (tr_A + D ** 0.5) / 2, (tr_A - D ** 0.5) / 2
  v_1, v_2 = [A[0][0] - lambda_2, A[1][0]], [A[0][0] - lambda_1, A[1][0]]

  return lambda_1, lambda_2, v_1, v_2
  
def eigenvectors_and_eigenvalues_second_type(A):
  tr_A, det_A = tr(A), det(A)
  D = determinant(1, -tr_A, det_A)

  lambda_1 = tr_A / 2

  return lambda_1

def eigenvectors_and_eigenvalues_third_type(A):
  tr_A, det_A = tr(A), det(A)
  D = determinant(1, -tr_A, det_A)

  Re_lambda = tr_A / 2
  Im_lambda = (-D) ** 0.5 / 2

  return Re_lambda, Im_lambda

def sum_1(a, b, n): # Re((a+ib)^n)
  coef = a ** n
  res = coef
  for i in range(1, n // 2 + 1):
    coef *= (-1) * b ** 2 / a ** 2 * (n - 2 * i + 2) * (n - 2 * i + 1) / (2 * i) / (2 * i - 1)
    res += coef
  return res

def sum_2(a, b, n): # Im((a+ib)^n)
  coef = n * a ** (n - 1) * b
  res = coef
  for i in range(2, (n + 1) // 2 + 1):
    coef *= (-1) * b ** 2 / a ** 2 * (n - 2 * i + 3) * (n - 2 * i + 2) / (2 * i - 1) / (2 * i - 2)
    res += coef
  return res

def find_c_first_type(lambda_1, lambda_2, v1, v2, a_0, b_0):
  c_1 = det([[a_0, v2[0]], [b_0, v2[1]]]) / det([[v1[0], v2[0]], [v1[1], v2[1]]])
  c_2 = det([[v1[0], a_0], [v2[1], b_0]]) / det([[v1[0], v2[0]], [v1[1], v2[1]]])
  return c_1, c_2

def find_c_second_type(lambda_1, a_0, b_0, x_0, x_1):
  c_1 = a_0
  a_1 = a_0 * x_0 + b_0 * x_1
  c_2 = a_1 - c_1 * lambda_1
  return c_1, c_2

def find_c_third_type(Re_lambda, Im_lambda, a_0, b_0, x_0, x_1, y_0, y_1):
  first_sum_1 = sum_1(Re_lambda, Im_lambda, 1)
  second_sum_1 = sum_2(Re_lambda, Im_lambda, 1)
  first_sum_2 = sum_1(Re_lambda, Im_lambda, 2)
  second_sum_2 = sum_2(Re_lambda, Im_lambda, 2)
  a_1, b_1 = a_0 * x_0 + b_0 * x_1, a_0 * y_0 + b_0 * y_1
  a_2 = a_1 * x_0 + b_1 * x_1
  c_1 = det([[a_1, second_sum_1], [a_2, second_sum_2]]) / det([[first_sum_1, second_sum_1], [first_sum_2, second_sum_2]])
  c_2 = det([[first_sum_1, a_1], [first_sum_2, a_2]]) / det([[first_sum_1, second_sum_1], [first_sum_2, second_sum_2]])
  return c_1, c_2

def _pow_with_mod(a, n): #Быстрое возведение в степень с взятием по модулю
  if n <= 0:
    return 1
  if n == 1:
    return a % INF
  if n % 2:
    return (a * (_pow_with_mod(a, n // 2) ** 2)) % INF
  return (_pow_with_mod(a, n // 2) ** 2) % INF

def pow_with_mod(a, n): #Учитывает отрицательные числа
  if a < 0 and n % 2:
    return - _pow_with_mod(abs(a), n)
  return _pow_with_mod(abs(a), n)

def solve_recur(n, a_0, b_0, x_0, x_1, y_0, y_1):
  A = [[x_0, x_1], [y_0, y_1]]
  type_rec = find_type_recur(A)

  if type_rec == 0:
    lambda_1, lambda_2, v_1, v_2 = eigenvectors_and_eigenvalues_first_type(A)
    c_1, c_2 = find_c_first_type(lambda_1, lambda_2, v_1, v_2, a_0, b_0)
    return (c_1 * pow_with_mod(lambda_1, n) * v_1[0] + c_2 * pow_with_mod(lambda_2, n) * v_2[0]) % INF

  elif type_rec == 1:
    lambda_1 = eigenvectors_and_eigenvalues_second_type(A)
    c_1, c_2 = find_c_second_type(lambda_1, a_0, b_0, x_0, x_1)
    return (c_1 * pow_with_mod(lambda_1, n) + c_2 * n * pow_with_mod(lambda_1, n - 1)) % INF

  else:
    Re_lambda, Im_lambda = eigenvectors_and_eigenvalues_third_type(A)
    c_1, c_2 = find_c_third_type(Re_lambda, Im_lambda, a_0, b_0, x_0, x_1, y_0, y_1)
    return (c_1 * sum_1(Re_lambda, Im_lambda, n) + c_2 * sum_2(Re_lambda, Im_lambda, n)) % INF
  
def main():
  t = int(input())

  for _ in range(t):
    n, a_0, b_0, x_0, x_1, y_0, y_1 = input().split()
    n = int(n)
    a_0, b_0, x_0, x_1, y_0, y_1 = float(a_0), float(b_0), float(x_0), float(x_1), float(y_0), float(y_1)
    print(solve_recur(n, a_0, b_0, x_0, x_1, y_0, y_1))

if __name__ == '__main__':
    main()
