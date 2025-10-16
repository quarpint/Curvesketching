from sage.all import *

@cached_function
def all_find_root(f, a=-1e11, b=1e11, nst_list = None):
    r"""
    Returns a nested list of all real roots of a function f using the built-in find_root().
    It kills your kernel for some functions like f = x^3.
    To extract an unnested list of the solutions you will need apply nested2list.
    
    INPUT:
    
    * "f" - sage.symbolic.Expression - function
    * "a" - sage.real_mpfr (default: -1e11) - lower bound
    * "b" - sage.real_mpfr (default: 1e11) - upper bound
    
    EXAMPLES:
    
        sage: all_find_root(x^2-1)
        [0.999999999999999, [-1.0]]
        sage: all_find_root(x^2+1)
        []
        sage: all_find_root(-3/17*x^5+x^2+x/2-1/8)
        [-0.6572541325361815, [1.9089620583045388, [0.1830545835378689]]]
    """
    if nst_list == None:
        nst_list = []
    if f(x=a) == 0:
        nst_list.append(a)
        nst_list.append(all_find_root(f, a+0.000001, b))
        return nst_list
    if f(x=b) == 0:
        nst_list.append(b)
        nst_list.append(all_find_root(f, a, b-0.000001))
        return nst_list
    try:
        nst = find_root(f, a, b)
        nst_list.append(nst)
        nst_list.append(all_find_root(f, a, nst-0.000001))
        nst_list.append(all_find_root(f, nst+0.000001, b))
        return [i for i in nst_list if i != []]
    except RuntimeError:
        return nst_list
    

def nested2list(f, nested):
    r"""
    Returns a rolled out, sorted version of a nested list. 
    Uses f.degree() to get all the zeros of polynomials, need a smarter iteration when implementing for all kinds of functions.
    
    INPUT:
    
    * "f" - sage.symbolic.Expression - function of variable x
    * "nested" - nasty nested list of length two or less
    
    EXAMPLES:
    
        sage: nested2list(-3/17*x^5+x^2+x/2-1/8, [-0.657, [1.908, [0.183]]])
        [-0.657000000000000, 0.183000000000000, 1.90800000000000]
    """
    roots = []
    if len(nested):
        for i in range(f.degree(x)):
            try:
                roots.append(nested[0])
                nested = nested[-1]
            except TypeError:
                break
    return sorted(roots)


def round3tuple(tup):
    r"""
    Takes a tuple and returns a rounded version - three digits after decimal points - of it. Used in CurveSketching.
    """
    return (round(tup[0], ndigits=3), round(tup[1], ndigits=3))


def text4plot(label, value, coord, col='blue', vert='bottom', rot = 0):
    r"""
    Returns a sage.plot.Graphics text element containing the label-value and the coordinates on the plot.
    Used in CurveSketching.
    
    INPUT:
    
    * "label" - str - text with a {} for each 'value' entry
    * "value" - many types - values to fill the label {} gaps
    * "coord" - tuple - coordinates of label-value on plot
    """
    return text(label.format(value), coord, fontsize = 'x-small', fontstyle= 'normal', horizontal_alignment='left',
                    vertical_alignment=vert, color=col, rotation = rot)