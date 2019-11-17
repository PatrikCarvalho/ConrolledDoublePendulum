
def A_lin(M1_val, M2_val, L1_val, L2_val, G_val, th_1_eql, th_1_d_eql, th_2_eql, th_2_d_eql):
    #import symengine
    import sympy as sym
    
    
    
    # Define variables and functions
    def Jacobian(v_str, c_str, f_list):
        c_str = sym.symbols(c_str)
        var = sym.symbols(v_str)
        f = sym.sympify(f_list)
        J = sym.zeros(len(f),len(var))
        for i, fi in enumerate(f):
            for j, s in enumerate(var):
                J[i,j] = sym.diff(fi, s)
        return J
    
    
    v_str = 'th_1 th_1_d th_2 th_2_d'
    c_str = 'M1 M2 L1 L2 G'
    
    f_1 = 'th_1_d' 
    f_2 = '(M2 * L1 * th_1_d * th_1_d * sin(th_2-th_1) * cos(th_2-th_1) + M2 * G * sin(th_2) * cos(th_2-th_1) + M2 * L2 * th_2_d * th_2_d * sin(th_2-th_1) - (M1+M2) * G * sin(th_1)) / ((M1+M2) * L1 - M2 * L1 * cos(th_2-th_1) * cos(th_2-th_1))'
    f_3 = 'th_2_d'
    f_4 = '(- M2 * L2 * th_2_d * th_2_d * sin(th_2-th_1) * cos(th_2-th_1) + (M1+M2) * G * sin(th_1) * cos(th_2-th_1) - (M1+M2) * L1 * th_1_d * th_1_d * sin(th_2-th_1) - (M1+M2) * G * sin(th_2)) / ((L2/L1) * ((M1+M2) * L1 - M2 * L1 * cos(th_2-th_1) * cos(th_2-th_1)))'
    f_list = [f_1, f_2, f_3, f_4]
    
    th_1, th_1_d, th_2, th_2_d, M1, M2, L1, L2, G = sym.symbols("th_1 th_1_d th_2 th_2_d M1 M2 L1 L2 G")
    
    J = Jacobian(v_str, c_str, f_list)
    #print(J.evalf(subs={M1:1, M2:1, L1:1, L2:1, G:9.81}))
    #print(J.simplify())
    #print(J)
    A = J.subs([(M1,M1_val), (M2,M2_val), (L1,L1_val), (L2,L2_val), (G,G_val), (th_1,th_1_eql), (th_1_d,th_1_d_eql), (th_2,th_2_eql), (th_2_d,th_2_d_eql)])
    return A.evalf(3)