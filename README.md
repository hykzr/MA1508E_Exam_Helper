# MA1505E Helper

A CLI app to solve common linear algebra questions in MA1508E Year 24/25, Sem2.

This srcipt solely depends on sympy. It has nothing to do with AI and is completely offline.

## Disclaimer

This script is provided for educational purposes only. The author is not responsible for any misuse of the script, including but not limited to academic dishonesty (e.g., cheating). The script may contain bugs, and the author is not liable for any damages caused. Users are responsible for verifying the results produced by the script.

## Features

- **Matrix Input Options**:
  - Manual entry (row-wise, space-separated values).
  eg:
  ```
  1 2 a
  3 4 x
  ```
  - Matlab Format (e.g., `[1,2;3,4]`)
  - SymPy Matrix string (e.g., `Matrix([[1, 2], [3, 4]])`).
  - Some LaTeX matrix or an array of vectors (e.g., `\begin{bmatrix} 1 & 2 \\ 3 & 4 \end{bmatrix}`).
  - if you need to input an array set you can input as a matrix where each column is one vector
- **Core Operations**:
  1. parse and format matrix (1+ will add omitted *)
  2. REF (this will display each step as well)
  3. RREF (input "3+" will display each step as well)
  4. Determinant
  5. Inverse
  6. Transpose
  7. Multiply Matrices
  8. Add Matrices
  9. Row Operations
  10. Solve A * X = B
  11. Solve X * A = B
  12. Find Vector Span
  13. Check if col Vector is in colspace
  14. Compare Subspaces
  15. Check if Vector Set is a Basis
  16. Substitude
  17. Null Space
  18. Check Orthogonal and Orthonormal Properties
  19. Orthogonal Projection
  20. Gram-Schmidt Orthogonalization
  21. Least Squares Solution
  22. Eigenvalues and Eigenvectors
  23. Diagnolization
  24. Equilibrium Vector
  25. differential systems
  26. intersection of subspaces
- **Output**:
  - pretty print (sympy.pprint) the matrix
  - print the result in matlab format
  - print the sympy string representation of the matrix
  - a set of vectors will converted to a matrix with each column is one vector

## Installation

1. **Prerequisites**:
   - Python 3.6 or higher.
   - SymPy library (`pip install sympy`).

2. clone the repository or simply download linear_algebra.py

## Usage

- Run the script and interact with the menu-driven interface:

```bash
python linear_algebra.py
```
- Run the code inside MatLab (there is something run with pretty print while running in a MatLab window)
```matlab
setenv('PYTHONIOENCODING', 'utf-8')
system("python linear_algebra.py")
```

## License

This project is licensed under the MIT License. See the [LICENSE](LICENSE) file for details.

## Acknowledgments

- Built with [SymPy](https://www.sympy.org/) for symbolic mathematics.
