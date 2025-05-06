import sympy as sp
import re
import os
from traceback import print_exc
global_symbols = {
    'Matrix': sp.Matrix,
    'sin':sp.sin, 'cos':sp.cos, 'tan':sp.tan, 'cot': sp.cot, 'sec': sp.sec, 'csc': sp.csc,
    'sqrt':sp.sqrt, 'pow': sp.Pow
}
def print_matrix(matrix, step, operation):
    print(f"Step {step}: {operation}")
    sp.pprint(matrix, use_unicode=True)
    print("\n" + "="*40 + "\n")
def preprocess_equation(eq):
    eq = re.sub(r'(?<=[a-zA-Z0-9)])(?=[a-zA-Z(])', '*', eq)
    return eq
def row_add_to_str(coe:sp.Expr, base_row, to_add_row):
    coe = sp.simplify(coe)
    if coe.is_Number:
        if coe == 1:
            return f"{base_row} + {to_add_row}"
        if coe == -1:
            return f"{base_row} - {to_add_row}"
        if coe > 0:
            return f"{base_row} + {coe} * {to_add_row}"
        if coe < 0:
            return f"{base_row} - {-coe} * {to_add_row}"
        return "[UNCHANGE]"
    elif coe.is_Add:
        return f"{base_row} + ({coe}) * {to_add_row}"
    else:
        return f"{base_row} + {coe} * {to_add_row}"
def row_multiply_to_str(coe:sp.Expr, base_row):
    coe = sp.simplify(coe)
    if coe == 1:
        return "[UNCHANGE]"
    elif coe.is_Add or coe.is_negative:
        return f"{base_row} * ({coe})"
    else:
        return f"{base_row} * {coe}"
def text_expr_is_constant(text_expr,variable_letters):
    for i in variable_letters:
        if i in text_expr:
            return False
    return True
def get_multiply_result_from_list(l,negative=False):
    if l==[]:
        if negative:
            return -1
        return 1
    return sp.simplify('*'.join(l))
def extract_coeff(text_expr:str,variable_letters):
    if not text_expr:
        return 0
    if text_expr[0]=='-':
        c,v=extract_coeff(text_expr[1:],variable_letters)
        return (-c,v)
    l=text_expr.split('*')
    coeff_expr=[]
    var_expr=[]
    for i in l:
        if text_expr_is_constant(i,variable_letters):
            coeff_expr.append(i)
        else:
            var_expr.append(i)
    return (get_multiply_result_from_list(coeff_expr),get_multiply_result_from_list(var_expr))
def parse_linear_system_to_expr_matrix(input_str:str, variable_letters=['x','y','z']):
    input_str = input_str.strip().replace('[', '').replace(']', '')
    if '\n' in input_str:
        lines = input_str.split('\n')
    else:
        lines = input_str.split(';')
    text_expr_matrix_left=[]
    text_expr_matrix_right=[]
    for line in lines:
        is_right=False
        text_expr_row_left=[]
        text_expr_row_right=[]
        curr_text_expr=''
        bracket_layer=0
        for ch in line:
            if ch==' ':
                continue
            if ch in ['(','[']:
                bracket_layer+=1
            elif ch in [')',']']:
                bracket_layer-=1
            if curr_text_expr!='':
                if ch in ['+','-','='] and bracket_layer==0:
                    curr_text_expr=extract_coeff(curr_text_expr.strip('+'),variable_letters)
                    if is_right:
                        text_expr_row_right.append(curr_text_expr)
                    else:
                        text_expr_row_left.append(curr_text_expr)
                    curr_text_expr=''
            if ch!='=' :
                curr_text_expr+=ch
            else:
                is_right=True
        if curr_text_expr!='':
            curr_text_expr=extract_coeff(curr_text_expr.strip('+'),variable_letters)
            if is_right:
                text_expr_row_right.append(curr_text_expr)
            else:
                text_expr_row_left.append(curr_text_expr)
        text_expr_matrix_left.append(text_expr_row_left)
        text_expr_matrix_right.append(text_expr_row_right)
    return (text_expr_matrix_left,text_expr_matrix_right)
def extract_all_vars(expr_matrix):
    vars=set()
    for row in expr_matrix:
        for _,var in row:
            if var and var!=1:
                vars.add(var)
    return vars
def parse_from_expr_matrix(text_expr_matrix_left, text_expr_matrix_right):
    matrix=[]
    vars=sorted(extract_all_vars(text_expr_matrix_left).union(extract_all_vars(text_expr_matrix_right)),key=lambda x: (len(str(x)),str(x)))
    print("variables: ",vars)
    for i in range(len(text_expr_matrix_left)):
        row=[0 for _ in range(len(vars)+1)]
        for coe,var in text_expr_matrix_left[i]:
            if var in vars:
                row[vars.index(var)]+=coe
            else:
                row[-1]-=coe
        for coe,var in text_expr_matrix_right[i]:
            if var in vars:
                row[vars.index(var)]-=coe
            else:
                row[-1]+=coe
        matrix.append(row)
    return matrix
def parse_latex_matrix(latex_str,add_mul):
    try:
        latex_str = latex_str.strip()
        # Match \begin{array}{spec} ... \end{array} or \left(\begin{array}{spec} ... \end{array}\right)
        matrix_pattern = r'\\begin\{array\}\{([clr| ]+)\}(.*?)\\end\{array\}'
        match = re.search(matrix_pattern, latex_str, re.DOTALL)
        if not match:
            # Try standard matrix environments (bmatrix, pmatrix, etc.)
            matrix_pattern_alt = r'\\begin\{[bpvV]?matrix\}(.*?)\\end\{[bpvV]?matrix\}'
            match = re.search(matrix_pattern_alt, latex_str, re.DOTALL)
            if not match:
                raise ValueError("Invalid LaTeX matrix format. Use \\begin{bmatrix} ... \\end{bmatrix} or \\begin{array}{ccc} ... \\end{array}.")
            content = match.group(1).strip()
        else:
            content = match.group(2).strip()
        rows = [row.strip() for row in content.split('\\\\') if row.strip()]
        if not rows:
            raise ValueError("Empty matrix content")
        
        matrix_entries = []
        cols = None
        for row in rows:
            if add_mul:
                row = preprocess_equation(row)
            elements = [elem.strip('{ }') for elem in row.split('&')]
            if cols is None:
                cols = len(elements)
            elif len(elements) != cols:
                raise ValueError("Inconsistent number of columns in matrix")
            parsed_row = [sp.sympify(elem, evaluate=True) for elem in elements]
            matrix_entries.append(parsed_row)
        
        return sp.Matrix(matrix_entries)
    except Exception as e:
        return None
def parse_latex_vector_set(latex_str,add_mul):
    try:
        latex_str = latex_str.strip()
        # Match individual vectors, e.g., \left(\begin{array}{l} ... \end{array}\right)
        vector_pattern = r'\\(left\()?\\begin\{array\}\{[clr]\}(.*?)\\end\{array\}\\?(right\))?'
        vectors = re.findall(vector_pattern, latex_str, re.DOTALL)
        if not vectors:
            raise ValueError("No valid vectors found in input. Use format like \\left(\\begin{array}{l} 1 \\\\ 2 \\end{array}\\right), ...")
        
        matrix_entries = []
        rows = None
        for _, content, _ in vectors:
            if add_mul:
                content = preprocess_equation(content)
            elements = [elem.strip() for elem in content.split('\\\\') if elem.strip()]
            if rows is None:
                rows = len(elements)
            elif len(elements) != rows:
                raise ValueError("Inconsistent number of rows in vector set")
            parsed_elements = [sp.sympify(elem, evaluate=True) for elem in elements]
            matrix_entries.append(parsed_elements)
        
        # Transpose to make each vector a column
        return list_to_col_matrix(matrix_entries)
    except Exception as e:
        print(f"Error parsing LaTeX vector set: {e}")
        return None
def parse_from_lateX(latex_str:str,add_mul:bool):
    try:
        return parse_latex_matrix(latex_str,add_mul)
    except:
        try:
            return parse_latex_vector_set(latex_str,add_mul)
        except:
            return None
def parse_matrix(input_str:str, add_mul:bool):
    latex_result=parse_from_lateX(input_str,add_mul)
    if latex_result:
        return latex_result
    if input_str.casefold().startswith('matrix'):
        return sp.simplify(input_str.replace('matrix','Matrix'))
    if input_str.startswith('['):
        try:
            return sp.simplify(f'Matrix({input_str})')
        except:
            pass
    input_str=input_str.strip('[').strip(']')
    if '=' in input_str:
        lhs,rhs=parse_linear_system_to_expr_matrix(input_str)
        return sp.Matrix(parse_from_expr_matrix(lhs,rhs))
    if '\n' in input_str:
        rows = input_str.split('\n')
    else:
        rows = input_str.split(';')
    matrix = []
    for row in rows:
        row = row.strip()
        if add_mul:
            row = preprocess_equation(row)
        if ',' in row:
            elements = row.split(',')
        elif ' ' in row:
            elements = row.split()
        else:
            try:
                matrix.append([sp.sympify(row)])
            except Exception as e:
                print("cannot parse:", elements, 'because:',repr(e))
            continue
        matrix.append([sp.sympify(e) for e in elements])
    return sp.Matrix(matrix)
def ref_with_steps(matrix):
    global last_output
    matrix = sp.Matrix(matrix)
    step = 0
    rows, cols = matrix.shape
    lead = 0
    
    for r in range(rows):
        if lead >= cols:
            break
        i = r
        while matrix[i, lead] == 0:
            i += 1
            if i == rows:
                i = r
                lead += 1
                if lead == cols:
                    break
        if lead >= cols:
            break
        if i != r:
            matrix.row_swap(i, r)
            step += 1
            print_matrix(matrix, step, f"R{r+1} <-> R{i+1}")
        
        for i in range(r + 1, rows):
            lv = matrix[i, lead]
            if lv != 0:
                coe = lv / matrix[r, lead]
                matrix.row_op(i, lambda v, j: v - coe * matrix[r, j])
                matrix = matrix.applyfunc(sp.simplify)
                step += 1
                print_matrix(matrix, step, row_add_to_str(-coe, f"R{i+1}", f"R{r+1}"))
        
        lead += 1
    last_output=matrix
    return matrix
def rref_with_steps(matrix):
    global last_output
    matrix = ref_with_steps(matrix)
    rows, cols = matrix.shape
    step = 0
    
    for r in range(rows - 1, -1, -1):
        lead_idx = None
        for j in range(cols):
            if matrix[r, j] != 0:
                lead_idx = j
                break
        if lead_idx is not None:
            lv = matrix[r, lead_idx]
            if lv!=1:
                matrix.row_op(r, lambda v, _: v / lv)
                matrix = matrix.applyfunc(sp.simplify)
                step += 1
                print_matrix(matrix, step, row_multiply_to_str(1/lv, f"R{r+1}"))
            for i in range(r):
                lv = matrix[i, lead_idx]
                if lv != 0:
                    matrix.row_op(i, lambda v, j: v - lv * matrix[r, j])
                    matrix = matrix.applyfunc(sp.simplify)
                    step += 1
                    print_matrix(matrix, step, row_add_to_str(-lv, f"R{i+1}", f"R{r+1}"))
    last_output=matrix
    return matrix

last_input=last_output=None
def input_matrix_helper(name = 'matrix', add_mul = False):
    print(f"Please input {name}, after inputing just press enter one more time\nyou may also input '_' for last result or '-' for last input")
    ipt=[]
    #read first line
    while 1:
        r = input().strip().replace('−','-').replace('','-')
        if not r:continue
        if r == '-':
            if last_input:
                return last_input
            print("ERROR: you have not inputted any matrix before")
        elif r in ['_', 'res', 'result']:
            if last_output:
                return last_output
            print("ERROR: there is no previous result")
        elif r[0].casefold() in 'i':
            try:
                nrows = int(r[1:])
            except ValueError:
                break
            if nrows == 0:
                print("Empty Matrix")
                return sp.Matrix([])
            return sp.Matrix(sp.Identity(nrows))
        elif r[0].casefold() == 'z':
            try:
                temp = r[1:].split()
                if len(temp) == 0: break
                if len(temp) == 1:
                    return sp.Matrix(sp.ZeroMatrix(int(temp[0]),int(temp[0])))
                elif len(temp) == 2:
                    return sp.Matrix(sp.ZeroMatrix(int(temp[0]),int(temp[1])))
                else:
                    break
            except ValueError:
                break
        else:
            break
    ipt.append(r)
    while 1:
        r=input().replace('−','-').replace('','-').strip()
        if not r:
            break
        ipt.append(r)
    return parse_matrix('\n'.join(ipt),add_mul)
def input_matrix(name='matrix', add_mul = False):
    global last_input
    last_input = input_matrix_helper(name, add_mul)
    return last_input
def print_matrix_matlab(matrix):
    # Prints the matrix in MATLAB syntax
    print("[", end="")
    for i in range(matrix.rows):
        print(", ".join(str(matrix[i, j]) for j in range(matrix.cols)), end="")
        if i < matrix.rows - 1:
            print("; ", end = '')
    print("]")
def clear():
    os.system('cls' if os.name == 'nt' else 'clear')

def determinant(matrix):
    try:
        return matrix.det()
    except sp.matrices.exceptions.NonSquareMatrixError:
        print("Error: Determinant can only be calculated for square matrices.")
        return None
def inverse(matrix):
    try:
        if matrix.det() == 0:
            print("Matrix is singular and cannot be inverted.")
            return None
        return matrix.inv()
    except sp.matrices.exceptions.NonSquareMatrixError:
        print("Error: Inverse can only be computed for square matrices.")
        return None
def transpose(matrix):
    return matrix.T
def multiply_matrices(mat1, mat2):
    try:
        return mat1 * mat2
    except ValueError:
        print("Error: Matrices cannot be multiplied due to incompatible dimensions.")
        return None
def add_matrices(mat1, mat2):
    try:
        return mat1 + mat2
    except ValueError:
        print("Error: Matrices must have the same dimensions to be added.")
        return None
def remove_double_spaces(s:str):
    t = []
    for i in s.split(' '):
        if i:
            t.append(i)
    return ' '.join(t)
def print_result(mat:sp.Matrix, change_last_output=True):
    global last_output
    if mat.rows==0:
        print("[EMPTY]")
        return
    sp.pprint(mat, use_unicode=True)
    print("MATLAB format: ", end = '')
    print_matrix_matlab(mat)
    print(f"python representation:", remove_double_spaces(repr(mat).replace('\n','')))
    if change_last_output:
        last_output=mat
def row_operations(matrix):
    global last_output
    print("Input row operations (e.g., R1+1/3*R2, R2/(a+1), R1<->R4). Type 'exit' to stop.")
    while True:
        row_vars = {f'R{i+1}': matrix.row(i) for i in range(matrix.rows)}
        print("Current Matrix:")
        sp.pprint(matrix, use_unicode=True)
        operation = input("Enter row operation: ").strip()
        if not operation:
            continue
        if operation.lower() in ['e','q','exit']:
            break
        operation = preprocess_equation(operation)
        print(operation)
        if '<->' in operation:
            try:
                row1, row2 = map(lambda x: int(x[1:].strip())-1, operation.split('<->'))
                print(f"swap(R{row1}, R{row2})")
                matrix.row_swap(row1, row2)
            except:
                print("Invalid swap operation format. Try again.")
        else:
            try:
                if '=' in operation:
                    target, expr = map(str.strip, operation.split('='))
                    row_idx = int(re.search(r'R(\d+)', target).group(1)) - 1
                    row_expr = sp.sympify(expr, locals=row_vars)
                    matrix[row_idx, :] = row_expr
                else:
                    row_expr = sp.sympify(operation, locals=row_vars)
                    row_idx = int(re.search(r'R(\d+)', operation).group(1)) - 1
                    matrix[row_idx, :] = row_expr
            except Exception as e:
                print(f"Invalid row operation format: {e}. Try again.")
    last_output=matrix
    return matrix
def solve_left_multiply(A, B):
    try:
        return A.pinv() * B
    except Exception as e:
        print(f"Error solving A * X = B: {e}")
        return None
def solve_right_multiply(A, B):
    try:
        return B * A.pinv()
    except Exception as e:
        print(f"Error solving X * A = B: {e}")
        return None
def vector_span(matrix):
    rref_matrix, pivots = matrix.rref()
    dimension = len(pivots)
    return dimension, matrix.cols
def is_vector_in_subspace(subspace_matrix, vector):
    print(vector)
    augmented_matrix = subspace_matrix.row_join(vector)
    rref_matrix, pivots = augmented_matrix.rref()
    return pivots[-1] < subspace_matrix.cols
def compare_subspaces(matrix1, matrix2):
    # Find the ranks of the subspaces
    rank1 = matrix1.rank()
    rank2 = matrix2.rank()
    # Check if matrix1 is contained in matrix2
    augmented1 = matrix2.row_join(matrix1)
    mat1_cont_mat2=augmented1.rank() == rank2
    augmented2 = matrix1.row_join(matrix2)
    mat2_cont_mat1=augmented2.rank() == rank1
    if mat1_cont_mat2 and mat2_cont_mat1:
        print("The two subspaces are the same")
    elif mat2_cont_mat1:
        print("The second subspace is contained within the first subspace.")
    elif mat1_cont_mat2:
        print("The second subspace is contained within the first subspace.")
    else:
        print("The two subspaces does not contain each other")
def is_basis(matrix):
    return matrix.rows == matrix.cols and matrix.rank() == matrix.rows
def null_space(matrix:sp.Matrix):
    try:
        null_vectors =  matrix.nullspace()
        return list_to_col_matrix(null_vectors)
    except Exception as e:
        print(f"Error computing null space: {e}")
        return None

def is_orthogonal(matrix):
    cols = matrix.cols
    for i in range(cols):
        for j in range(i + 1, cols):
            if matrix[:, i].dot(matrix[:, j]) != 0:
                return False
    return True

def is_orthonormal(matrix):
    if not is_orthogonal(matrix):
        return False
    for i in range(matrix.cols):
        if matrix[:, i].norm() != 1:
            return False
    return True

def orthogonal_projection(vector, subspace_matrix):
    try:
        U = subspace_matrix
        proj = U * (U.T * U).inv() * U.T * vector
        return proj
    except Exception as e:
        print(f"Error computing orthogonal projection: {e}")
        return None
def list_to_col_matrix(vectors: list[sp.Matrix]) -> sp.Matrix:
    result = []
    for i in vectors:
        if i.cols == 1:
            i = i.T #convert a col vector to a row vector
        if i.rows != 1:
            print(f"Cannot Convert {i} to a colomn as it is not a vector")
        result.append(i)
    return sp.Matrix(result).T
def gram_schmidt(matrix: sp.Matrix, normalize=True):
    try:
        orthogonal_vectors = []
        for i in range(matrix.cols):
            vec = matrix[:, i]
            for prev in orthogonal_vectors:
                proj = (vec.dot(prev) / prev.dot(prev)) * prev
                vec -= proj
            orthogonal_vectors.append(vec)
        if normalize:
            for v in orthogonal_vectors:
                if v.norm()==0:
                    print("Cannot Normalize Matrix as is it not linearly independent")
                    return (list_to_col_matrix(orthogonal_vectors), None)
            orthonormal_vectors = [v / v.norm() for v in orthogonal_vectors]
            return (list_to_col_matrix(orthogonal_vectors), list_to_col_matrix(orthonormal_vectors))
        return (list_to_col_matrix(orthogonal_vectors), None)
    except Exception as e:
        print(f"Error in Gram-Schmidt process: {e}")
        return None
def extract_coefficients(solution: sp.Matrix, params: sp.Matrix):
    """
    Extracts coefficients of parameters and constant terms from a SymPy least_squares solution.

    Args:
        solution (Matrix): SymPy Matrix containing the solution expressions.
        params (Matrix): SymPy Matrix containing the parameter symbols (e.g., tau0, tau1).

    Returns:
        dict: A dictionary where:
              - keys are parameter names (str) and 'constant',
              - values are coefficient vectors (as SymPy Matrices).
    """
    param_list = list(params)
    coeff_dict = {param: [] for param in param_list}
    coeff_dict['constant'] = []

    for expr in solution:
        # For each parameter, extract the coefficient
        for param in param_list:
            coeff = expr.coeff(param)
            coeff_dict[param].append(coeff)
        # Extract constant term (by removing all param parts)
        constant = expr
        for param in param_list:
            constant -= expr.coeff(param) * param
        coeff_dict['constant'].append(constant)

    # Convert each list into a SymPy Matrix
    for key in coeff_dict:
        coeff_dict[key] = sp.Matrix(coeff_dict[key])

    return coeff_dict
def least_squares_solution(A:sp.Matrix, b:sp.Matrix):
    b = A.T*b
    A = A.T*A
    return A.gauss_jordan_solve(b)
def eigenvalues_eigenvectors(matrix:sp.Matrix):
    try:
        if matrix.rows != matrix.cols:
            raise ValueError("Matrix must be square for eigenvalue/eigenvector computation")
        eigenvals = matrix.eigenvals()  # Dictionary of eigenvalue: multiplicity
        eigenvects = matrix.eigenvects()  # List of (eigenvalue, multiplicity, [eigenvectors])
        return eigenvals, eigenvects
    except Exception as e:
        print(f"Error computing eigenvalues and eigenvectors: {e}")
        print_exc()
        return None, None
def process_eigen():
    mat = input_matrix()
    if mat is not None:
        if mat.rows != mat.cols:
            print("Error: Matrix must be square for eigenvalue/eigenvector computation.")
        else:
            eigenvals, eigenvects = eigenvalues_eigenvectors(mat)
            if eigenvals is not None and eigenvects is not None:
                print("Eigenvalue [algebraic multiplicity]: {", end = ' ')
                t = []
                for val, mult in eigenvals.items():
                    t.append(f"{val}[{mult}]")
                print(', '.join(t),'}\nEigenvectors:')
                for val, alg_mult, vecs in eigenvects:
                    geom_mult = len(vecs)  # Geometric multiplicity is the number of eigenvectors
                    print(f"For eigenvalue {val}:")
                    print(f"Algebraic Multiplicity: {alg_mult}, Geometric Multiplicity: {geom_mult}")
                    print_result(list_to_col_matrix(vecs))
def check_diagonalizable(A:sp.Matrix):
    if not A.rows == A.cols:
        print("Not a square matrix so it is not diagonalizable")
    eigen_data = A.eigenvects()
    P_cols = []
    for val, mult, basis in eigen_data:
        if len(basis) < mult:
            print("Matrix is NOT diagonalizable.")
            return
        P_cols.extend(basis)
    P = sp.Matrix.hstack(*P_cols)
    D = sp.Matrix.zeros(A.shape[0])
    # Fill D with eigenvalues matching P's columns
    col_idx = 0
    for val, mult, basis in eigen_data:
        for _ in basis:
            D[col_idx, col_idx] = val
            col_idx += 1
    print("Matrix is diagonalizable.")
    print("P =")
    print_result(P)
    print("D =")
    print_result(D)
    
    # Check if matrix is symmetric (then automatically orthogonally diagonalizable)
    if A == A.T:
        print("Matrix is symmetric. It is orthogonally diagonalizable.")
    else:
        print("Matrix is NOT symmetric. It is NOT orthogonally diagonalizable")
        return

    # Try Gram-Schmidt within eigen spaces
    new_cols = []
    for val, mult, basis in eigen_data:
        B = sp.Matrix.hstack(*basis)
        _, B_orth = gram_schmidt(B)
        for i in range(B_orth.cols):
            new_cols.append(B_orth.col(i))

    P_orth = list_to_col_matrix(new_cols)
    print("After Gram-Schmidt orthogonalization:")
    print_result(P_orth)
def equilibrium_vector(S):
    try:
        if S.rows != S.cols:
            raise ValueError("Matrix S must be square")
        for i in range(S.rows):
            for j in range(S.cols):
                if S[i, j] < 0:
                    raise ValueError("Matrix S must have non-negative entries")
        is_stochastic = True
        S_stoch = S.copy()
        for i in range(S.cols):
            col_sum = sum(S[j, i] for j in range(S.rows))
            if col_sum!=1:
                is_stochastic = False
            if col_sum == 0:
                raise ValueError(f"Col {i+1} sums to zero, cannot normalize")
            for j in range(S.rows):
                S_stoch[j, i] = S[j, i] / col_sum
        
        n = S_stoch.rows
        S_power = S_stoch.copy()
        is_regular = False
        for k in range(1, n + 1):
            all_positive = all(S_power[i, j] > 0 for i in range(n) for j in range(n))
            if all_positive:
                is_regular = True
                break
            S_power = S_power * S_stoch
        I = sp.eye(n)
        A = S_stoch - I
        constraint = sp.Matrix([[1]*n])  # Row of ones
        A_aug = A.col_join(constraint)
        b = sp.Matrix([0]*(n) + [1])  # Zeroes except last entry = 1
        # Solve using pseudo-inverse for robustness
        if b.cols!=1 and b.rows==1:
            b=b.T
        x = A_aug.pinv() * b
        # Normalize to ensure sum = 1 (in case of numerical errors)
        x_sum = sum(x[i] for i in range(n))
        x = x / x_sum
        return S_stoch, is_stochastic, is_regular, x
    except Exception as e:
        print(f"Error computing equilibrium vector: {e}")
        return None, None, None, None
def  process_equilibrium_vector(S: sp.Matrix):
    if S.rows != S.cols:
        print("Error: Matrix S must be square.")
    else:
        S_stoch, is_stochastic, is_regular, x = equilibrium_vector(S)
        if S_stoch is not None:
            if not is_stochastic:
                print("Matrix S was not stochastic. Converted to stochastic matrix:")
                print_result(S_stoch, change_last_output=False)
            if not is_regular:
                print("Matrix is not regular, NO equilibrium vector")
            else:
                print_result(x)
def real_exp(eigenvalue:sp.Expr, t:sp.Symbol):
    real,image = eigenvalue.as_real_imag()
    return (sp.cos(image)+sp.I*sp.sin(image))*sp.exp(real*t)
# def solve_ode_system(A):
#     t = sp.Symbol('t')
#     n = A.shape[0]
#     # Compute eigenvalues and eigenvectors
#     P, J = A.jordan_form()
#     sp.pprint(P)
#     sp.pprint(J)
#     # Initialize solution vector
#     sol = sp.zeros(n, 1)
    
#     # Process each Jordan block
#     block_start = 0
#     while block_start < n:
#         eigenvalue = J[block_start, block_start]
#         # Determine size of Jordan block
#         k = 0
#         while (block_start + k + 1 < n and 
#                J[block_start + k, block_start + k + 1] == 1):
#             k += 1
#         k = k + 1  # Size of the block
        
#         # Compute the term for this block
#         zero_vec = sp.zeros(n, 1)
#         coe_mat = [(t**m / sp.factorial(m)) * P[:, block_start + k - 1 - m] for m in range(k)]
#         print(coe_mat)
#         term = real_exp(eigenvalue, t) * sum(
#             [(t**m / sp.factorial(m)) * P[:, block_start + k - 1 - m] for m in range(k)], zero_vec
#         )
#         sol += term
#         block_start += k
    
#     return sol
def get_particular_solution(N, rhs):
    """
    Solve N*x = rhs and return one particular solution.
    If free parameters exist, set them to 0.
    """
    # sp.linsolve returns a set of solutions in parametric form.
    sol_set = sp.linsolve((N, rhs))
    if not sol_set:
        # No solutions exist.
        return None
    # Pick the first solution (a tuple of expressions).
    sol = list(sol_set)[0]
    # Force a column vector.
    x_part = sp.Matrix(sol)
    to_subs = {symbol: 0 for symbol in x_part.free_symbols}
    x_part=x_part.subs(to_subs)
    return x_part

def solve_linear_system(A):
    """Solve the system y' = A y and return the general solution."""
    if A.rows != A.cols:
        raise ValueError("Matrix A must be square")
    
    t = sp.symbols('t')
    n = A.rows
    solution = []
    
    # Get eigenvalues and eigenvectors
    eigenvects = A.eigenvects()
    
    for eigenvalue, multiplicity, eigenvectors in eigenvects:
        if eigenvalue.is_real:
            # Real eigenvalues
            for vec in eigenvectors:
                solution.append(sp.exp(eigenvalue * t) * vec)
            if multiplicity > len(eigenvectors):
                vec = eigenvectors[0]
                vec_chain = []
                mat: sp.Matrix = (A - eigenvalue * sp.eye(n))
                for k in range(1, multiplicity - len(eigenvectors) + 1):
                    print('General Vector for ')
                    sp.pprint(vec)
                    vec = get_particular_solution(mat, vec)
                    sp.pprint(vec)
                    vec_chain.append(vec)
                    col = t**k*eigenvectors[0]/sp.factorial(k)
                    for i in vec_chain:
                        k-=1
                        col+=(t**k*i/sp.factorial(k))
                    solution.append(sp.exp(eigenvalue * t)*col)
        else:
            # Complex eigenvalues (λ = a + bi)
            a = sp.re(eigenvalue)
            b = sp.im(eigenvalue)
            if b<0:
                continue
            for vec in eigenvectors:
                u = sp.re(vec)  # Real part of eigenvector
                v = sp.im(vec)  # Imaginary part of eigenvector
                # Real solution from complex conjugate pair
                #print(a,b,u,v)
                solution.append(sp.exp(a * t) * (sp.cos(b * t) * u - sp.sin(b * t) * v))
                solution.append(sp.exp(a * t) * (sp.sin(b * t) * u + sp.cos(b * t) * v))
    
    # Form the general solution with arbitrary constants
    c = sp.symbols(f'c1:{len(solution)+1}')
    #general_solution = sum(c[i] * sol for i, sol in enumerate(solution))
    return list_to_col_matrix(solution)
def print_eq(lhs, rhs):
    l = sp.pretty(lhs)
    r = sp.pretty(rhs)
    if '\n' in l:
        print(l)
    else:
        print(l,'= ', end = '')
    if '\n' in r and not '\n' in l:
        print()
    print(r)
def process_ode(mat:sp.Matrix):
    res = solve_linear_system(mat)
    print("result set")
    print_result(res)
    for i in range(res.rows):
        y_symbol = sp.Symbol(f'y{i}')
        y = 0
        for j in range(res.cols):
            c_symbol = sp.Symbol(f'C{i}')
            y += c_symbol*res[i,j]
        print_eq(y_symbol, y)
        print('When t = 0')
        print_eq(y_symbol, y.subs('t',0))
def subspace_intersection(A, B):
    try:
        if A.rows != B.rows:
            raise ValueError("Matrices A and B must have the same number of rows")
        # Form the matrix [A | -B]
        neg_B = -B
        M = A.row_join(neg_B)
        # Compute the null space of [A | -B]
        null_vectors = M.nullspace()
        # Extract vectors in the intersection
        intersection_basis = []
        m = A.cols
        for vec in null_vectors:
            u = sp.Matrix(vec[:m])  # Coefficients for A
            x = A * u    # Vector in col(A), also equals B * v where v = vec[m:]
            # Only include non-zero vectors
            if not x.is_zero_matrix:
                intersection_basis.append(x)
        # Ensure linear independence (SymPy's nullspace typically returns a basis)
        return list_to_col_matrix(intersection_basis)
    except Exception as e:
        print(f"Error computing subspace intersection: {e}")
        print_exc()
        return None
def process_input():
    print("""Linear Algebra
| 0. exit                                                    |
| 1. parse and format matrix (1+ will add omitted *)         |
| 2. REF (this will display each step as well)               |
| 3. RREF (input "3+" will display each step as well)        |
| 4. Determinant                                             |
| 5. Inverse                                                 |
| 6. Transpose                                               |
| 7. Multiply Matrices                                       |
| 8. Add Matrices                                            |
| 9. Row Operations                                          |
|10. Solve A * X = B                                         |
|11. Solve X * A = B                                         |
|12. Find Vector Span                                        |
|13. Check if col Vector is in colspace                      |
|14. Compare Subspaces                                       |
|15. Check if Vector Set is a Basis                          |
|16. Substitude                                              |
|17. Null Space                                              |
|18. Check Orthogonal and Orthonormal Properties             |
|19. Orthogonal Projection                                   |
|20. Gram-Schmidt Orthogonalization                          |
|21. Least Squares Solution                                  |
|22. Eigenvalues and Eigenvectors                            |
|23. Diagnolization                                          |
|24. Equilibrium Vector                                      |
|25. differential systems                                    |
|26. intersection of subspaces                               |
+------------------------------------------------------------+""")
    while True:
        option = input("please choose one number: ").casefold().strip()
        if option:
            break
    if option in ['0', 'e', 'exit']:
        return False
    elif option in ['1','1+','parse','parse+']:
        if option[-1] == '+':
            mat = input_matrix(add_mul = True)
        else:
            mat = input_matrix()
        print_result(mat)
    elif option == '2' or option =='ref':
        mat = input_matrix()
        ref_with_steps(mat)
    elif option == '3' or option == 'rref':
        mat = input_matrix()
        print("RREF result")
        print_result(mat.rref()[0])
    elif option == '3+' or option == 'rref+':
        mat = input_matrix()
        rref_with_steps(mat)
    elif option == '4' or option == 'det':
        mat = input_matrix()
        det=determinant(mat)
        if det:
            print("Determinant:", det)
    elif option == '5' or option == 'inv':
        mat = input_matrix()
        inv_mat = inverse(mat)
        if inv_mat is not None:
            print("Inverse Matrix:")
            print_result(inv_mat)
    elif option == '6' or option == 'transpose':
        mat = input_matrix()
        print("Transpose Matrix:")
        print_result(transpose(mat))
    elif option == '7' or option == 'mul' or option == '*':
        mat1 = input_matrix("matrix A")
        mat2 = input_matrix("matrix B")
        result = multiply_matrices(mat1, mat2)
        if mat1==mat2:
            print("A^2 = ")
            print_result(result)
        else:
            if result is not None:
                print("A * B = ")
                print_result(result)
            result = multiply_matrices(mat2, mat1)
            if result is not None:
                print("B * A = ")
                print_result(result,False)
    elif option == '8' or option == 'add' or option == '+':
        mat1 = input_matrix("first matrix")
        mat2 = input_matrix("second matrix")
        result = add_matrices(mat1, mat2)
        if result is not None:
            print("Matrix Addition Result:")
            print_result(result)
    elif option == '9' or option == 'row':
        mat = input_matrix()
        row_operations(mat)
    elif option == '10' or option == 'sright':
        matA = input_matrix("A")
        matB = input_matrix("B")
        result = solve_left_multiply(matA, matB)
        if result is not None:
            print("Solution for A * X = B:")
            print_result(result)
    elif option == '11' or option == 'sleft':
        matA = input_matrix("A")
        matB = input_matrix("B")
        result = solve_right_multiply(matA, matB)
        if result is not None:
            print("Solution for X * A = B:")
            print_result(result)
    elif option == '12' or option == 'span':
        mat = input_matrix()
        dimension, space = vector_span(mat)
        if dimension!=space:
            print(f"The vector span is a {dimension}-dimensional subspace of R^{space}.")
        else:
            print(f"The vectors span the whole R^{dimension}")
    elif option == '13' or option == 'vec':
        subspace = input_matrix("subspace matrix")
        vec = input_matrix("vector")
        for i in range(vec.cols):
            result = is_vector_in_subspace(subspace, vec.col(i))
            print(f"The vector {vec.col(i)} is in the subspace." if result else f"The vector {vec.col(i)} is NOT in the subspace.")
    elif option == '14' or option=='cmpspan':
        mat1 = input_matrix("first subspace matrix")
        mat2 = input_matrix("second subspace matrix")
        compare_subspaces(mat1, mat2)
    elif option == '15' or option == 'isbasis':
        mat = input_matrix()
        result = is_basis(mat)
        print("The vector set forms a basis." if result else "The vector set does NOT form a basis.")
    elif option == '16' or option == 'subs':
        mat = input_matrix()
        subexpr=preprocess_equation(input("please input the substitution (eg. a=1, x=2, x=y+2): ").strip()).split("=")
        if len(subexpr)!=2:
            print("Invalid substitution equation")
        try:
            print_result(mat.subs(subexpr[0].strip(),subexpr[1].strip()))
        except Exception as e:
            print(f"ERROR while substituting {subexpr} :",e)
    elif option == '17' or option == 'null':
        mat = input_matrix()
        if mat is not None:
            null_vecs = null_space(mat)
            if null_vecs is not None:
                print("Null Space Basis Vectors (one col one vector):")
                print_result(null_vecs)
    elif option == '18' or option == 'checkortho':
        mat = input_matrix()
        if mat is not None:
            ortho = is_orthogonal(mat)
            orthonorm = is_orthonormal(mat)
            if orthonorm:
                print("The set is orthonormal (and thus orthogonal).")
            elif ortho:
                print("The set is orthogonal but not orthonormal.")
            else:
                print("The set is neither orthogonal nor orthonormal.")
    elif option == '19' or option == 'proj':
        subspace_mat = input_matrix("subspace matrix")
        vec = input_matrix("vectors to project")
        if subspace_mat is not None and vec is not None:
            if vec.rows != subspace_mat.rows:
                print("Error: Vectosr must have the same number of rows as the subspace matrix.")
            else:
                for i in range(vec.cols):
                    proj = orthogonal_projection(vec.col(i), subspace_mat)
                    if proj is not None:
                        print("Orthogonal Projection:")
                        print_result(proj)
    elif option == '20' or option == 'gs':
        mat = input_matrix()
        if mat is not None:
            orthognal_mat, orthonamal_mat = gram_schmidt(mat)
            if orthognal_mat:
                print("Orthognal Basis after Gram-Schmidt:")
                print_result(orthognal_mat)
            if orthonamal_mat:
                print("Orthonormal Basis after Gram-Schmidt:")
                print_result(orthonamal_mat)
    elif option == '21' or option == 'ls':
        A = input_matrix("coefficient matrix A")
        b = input_matrix("vector b")
        if A is not None and b is not None:
            if b.cols != 1 or b.rows != A.rows:
                print("Error: Vector b must be a column vector with the same number of rows as matrix A.")
            else:
                solution, params = least_squares_solution(A, b)
                if solution is not None:
                    print("Least Squares Solution:")
                    print_result(solution)
                if params:
                    print('Parameters Coefficients')
                    for key, vec in extract_coefficients(solution, params).items():
                        if type(key) == str:
                            print(f"{key}: ", end='')
                        else:
                            print(f'{sp.pretty(key)}: ', end = '')
                        sp.pprint(vec.T, use_unicode=True)
                        
    elif option == '22' or option == 'eigen':
        process_eigen()
    elif option in ['23', 'dia']:
        check_diagonalizable(input_matrix())
    elif option in ['24', 'eqvec']:
        process_equilibrium_vector(input_matrix("transition matrix S"))
    elif option in ['25', 'ode']:
        process_ode(input_matrix("Coefficient Matrix A for y' = Ay"))
    elif option in ['26', 'inter']:
        A = input_matrix("matrix A (first subspace)")
        B = input_matrix("matrix B (second subspace)")
        if A is not None and B is not None:
            if A.rows != B.rows:
                print("Error: Matrices A and B must have the same number of rows.")
            else:
                basis = subspace_intersection(A, B)
                if basis is not None:
                    if not basis:
                        print("The intersection of the subspaces is {0}.")
                    else:
                        print("Basis for the intersection of the subspaces:")
                        print_result(basis)
    else:
        print("???")

                
    input("Press enter to continue")
    return True
def main():
    while True:
        try:
            if not process_input():
                return
        except KeyboardInterrupt:
            print("ctrl+C detected, exiting")
            break
        except Exception as e:
            print("ERROR: ", repr(e))
            print_exc()
            input("Press enter to continue")
if __name__ == '__main__':
    main()