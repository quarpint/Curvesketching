import random
from sage.all import *
from .workarounds import *


class CurveSketching:
    r"""
    A class for curve sketching.
    """
    def __init__(self, f, a=-1e11, b=1e11, x0=0, method='solve'):
        """
        All the necessary attributes for CurveSketching (on an interval if desired). I recommend using 'method' roots for large polynomials.
        One could add more attributes e.g. zeros, extrema, integral, tangent resp. inherit from classes for much more fun.
        
        INPUT:
        
        * "f" - sage.symbolic.Expression - function (of variable x needed for most features)
        * "a" - sage.Integer/real_mpfr/... (default: -1e11) - lower bound of interval
        * "b" - sage.Integer/real_mpfr/... (default: 1e11) - upper bound of interval
        * "x0" - sage.Integer/real_mpfr/... (default: 0) - xcoordinate to calculate the tangent
        * "method" - 'solve' or 'roots' (default: 'solve') - get the zeros (incl. extrema) with the solve() or roots() method respectively
        
        EXAMPLES:
            
            CurveSketching(x^123+3*x^312+7*x+1, 0, 2, x0=2, method='roots')
            CurveSketching(x^2 - 1)
            %prun -s cumulative print(CurveSketching(3*x^312 + x^123 + 7*x + 1, method='roots'))
            %prun -s cumulative print(CurveSketching(3*x^312 + x^123 + 7*x + 1, method='solve'))
        """
        var('x')
        self.f = f
        self.a = a
        self.b = b
        self.x0 = x0
        self.fx = diff(f, x)
        self.method = method
    

    def zeros_solve(self, fx = None):
        r"""
        Returns a sorted list with the xcoordinates of all real roots using the built-in solve().
        Takes a optional function input 'fx' to calculate the roots of. 
        
        EXAMPLES:
            
            sage: CurveSketching(x^2-1, 0, 2).zeros_solve()
            [1.00000000000000]
            sage: CurveSketching(-3/17*x^5+x^2+x/2-1/8).zeros_solve()
            [-0.657254138266796, 0.183054585554126, 1.90896209705987]
            sage: CurveSketching(x^2+1).zeros_solve()
            []
        """
        if fx == None:
            roots = solve(self.f, x, to_poly_solve = 'force')
        else:
            roots = solve(fx, x, to_poly_solve = 'force')
        return sorted([n(i.rhs()) for i in roots if i.rhs().imag()==0 and i.rhs()>=self.a and i.rhs()<=self.b])

    
    def zeros_roots(self, fx = None):
        r"""
        Returns a sorted list with the xcoordinates of all real roots using the built-in roots() (specific for polynomials).
        Takes a optional function input 'fx' to calculate the roots of. 
        
        EXAMPLES:
            
            sage: CurveSketching(x^2-1, 0, 2).zeros_roots()
            [1.00000000000000]
            sage: CurveSketching(-3/17*x^5+x^2+x/2-1/8).zeros_roots()
            [-0.657254132536182, 0.183054583537869, 1.90896205830454]
            sage: CurveSketching(x^2-1).zeros_roots(x^2+1)
            []
        """
        if fx == None:
            roots = self.f.roots(ring = RR)
        else:
            roots = fx.roots(ring=RR)
        return [n(i[0]) for i in roots if i[0]>=self.a and i[0]<=self.b]

    
    def zeros_find_root(self):
        r"""
        Returns a sorted list with the xcoordinates of all real roots using all_find_root (via built-in find_root - num approx).
        
        EXAMPLES:
            
            sage: CurveSketching(x^2-1, 0, 2).zeros_find_root()
            [1.0]
            sage: CurveSketching(-3/17*x^5+x^2+x/2-1/8).zeros_find_root()
            [-0.6572541325361815, 0.1830545835378689, 1.9089620583045388]
            sage: CurveSketching(x^2+1).zeros_find_root()
            []
        """
        nested = all_find_root(self.f, self.a, self.b)
        return nested2list(self.f, nested)

    
    def zeros_newton(self, epsilon = 1/1e11):
        r"""
        Returns a sorted list with the xcoordinates of all real roots using a sloppy self implemented newton method (num. approx).
        Takes an optional argument 'epsilon' (default: 1/1e11) as stopping precision.
        
        EXAMPLES:
            
            sage: CurveSketching(x^2-1, 0, 2).zeros_newton()
            [1.0000000000000000000]
            sage: CurveSketching(-3/17*x^5+x^2+x/2-1/8).zeros_newton()
            [-0.65725413253618114595, 0.18305458353786863146, 1.9089620583045387383]
            sage: CurveSketching(x^2+1).zeros_newton()
            []
        """
        roots = []
        func = self.f
        for i in range(self.f.degree(x)):
            funcx = diff(func,x)
            newton = x - (func/funcx)(x=x)
            
            random.seed(69)
            xn = random.uniform(self.a, self.b) #starting point, careful if it is an extremum
            while funcx(x=xn) == 0:
                xn = random.uniform(self.a, self.b) + abs(epsilon)
            
            stop = [Infinity, xn]
            while abs(stop[-1] - stop[-2]) > abs(epsilon): #till difference between consecutive elements gets small
                xn = N(newton(x=xn), digits=20)
                stop.append(xn)
                
                if len(stop) == 1000: #we assume we found every root or there is none
                    return sorted([i for i in roots if i>=self.a and i<=self.b])
            
            roots.append(xn)
            g = x - xn
            func = func.maxima_methods().divide(g)[0]
        
        return sorted([i for i in roots if i>=self.a and i<=self.b])
 

    def derivative(self, order=1):
        r"""
        Returns derivative. Or derivative of a chosen 'order' (default: 1).
        """
        if order == 1:
            return self.fx
        else:
            return diff(self.f, x, order)
 

    def extrema(self):
        r"""
        Returns a tuple containing three lists of the minima, maxima and inflection points respectively.
        Uses the solve() method to calculate the zeros.
        
        EXAMPLES:
            
            sage: CurveSketching(x^2-1, 0.1, 2).extrema()
            ([], [], [])
            sage: CurveSketching((x+17)*(x-3)*(x-1/8)^3).extrema()
            ([2.30095232337017], [-13.4509523233702], [0.125000000000000])
            sage: CurveSketching(x^2+1).extrema()
            ([0.000000000000000], [], [])
        """
        if self.method == 'solve':
            candidates = self.zeros_solve(fx = self.fx)
        else:
            candidates = self.zeros_roots(fx = self.fx)
        
        fxx = diff(self.fx, x)
        minima = []
        maxima = []
        inflection = []
        
        if fxx.degree(x) != 0:
            for e in candidates:
                if fxx(x=e) < 0:
                    maxima.append(e)
                if fxx(x=e) > 0:
                    minima.append(e)
                if fxx(x=e) == 0:
                    fxxx = diff(fxx, x)
                    if fxxx(x=e) != 0:
                        inflection.append(e)
        else:
            for e in candidates:
                if fxx < 0:
                    maxima.append(e)
                if fxx > 0:
                    minima.append(e)
                if fxx == 0:
                    print('There is a non-identifiable extremum.')
        return minima, maxima, inflection

    
    def integration(self):
        r"""
        Solves the integral. If no interval is specified it returns the indefinite integral (without constant C).
        
        EXAMPLES:
            
            sage: CurveSketching(x^2-1, 0, 2).integration()
            2/3
            sage: CurveSketching((x+17)*(x-3)*(x-1/8)^3).integration()
            1/6*x^6 + 109/40*x^5 - 3597/256*x^4 + 10127/1536*x^3 - 619/512*x^2 + 51/512*x
            sage: CurveSketching(sin(x)*x*(3+cos(x))).integration()
            -1/4*x*cos(2*x) - 3*x*cos(x) + 1/8*sin(2*x) + 3*sin(x)
        """
        try:
            if abs(self.a) == self.b == 1e11:
                return integral(self.f, x)
            return integral(self.f, (x, self.a, self.b))
        except ValueError:
            return None

        
    def tangent(self, x0_new=None):
        r"""
        Returns the equation of the tangent at a specified point 'x0' (default: 0).
        
        EXAMPLES:
        
            sage: CurveSketching(x^2-1, 0, 2, x0=2).tangent()
            4*x - 5
            sage: CurveSketching(x^2-1).tangent()
            -1
            sage: CurveSketching(x^2-1).tangent(2)
            4*x - 5
        """
        if x0_new != None:
            self.x0 = x0_new
        m = diff(self.f, x)(x=self.x0)
        g = m*(x-self.x0) + self.f(x=self.x0)
        return g
 

    def vertical_asymptotes(self):
        r"""
        Returns the xcoordinates of the vertical asymptotes.
        """
        ans = solve(1/self.f, x)
        return sorted([n(i.rhs()) for i in ans])

    
    def horizontal_asymptotes(self, a=-Infinity, b=Infinity):
        r"""
        Returns the ycoordinates of the vertical asymptotes. Takes optionally 'a' and 'b' if you want the interval limits.
        """
        return limit(self.f, x=a), limit(self.f, x=b)
  

    def plot(self, xmin=-10, xmax=10, ymin=-10, ymax=10, tan = False):
        r"""
        Returns a plot with the labeled function, zeros (solve() method), extreme values resp. inflection points, tangent and asymptotes.
        If you have difficulties displaying the graphics try on a smaller interval.
        
        INPUT:
        
        * "xmin", "xmax" - sage.Integer/real_mpfr/... (default: -10 , 10) - x coordinate range of plot
        * "ymin", "ymax" - sage.Integer/real_mpfr/... (default: -10 , 10) - y coordinate range of plot
        * "tan" - True or False (default: False) - set to True to add the tangent to the plot (careful it can change your plotting range)
        
        EXAMPLES:
        
            CurveSketching(x/(x^2-1)+1, 0, 2, x0=2).plot(tan=True)
            CurveSketching(-3/17*x^5+x^2+x/2-1/8, -12, 20, x0=2).plot(-5,5,-100,100,tan=True)
            CurveSketching(sin(x)^2, 0, 4*pi, x0=1/sqrt(2)).plot(0,4*pi, tan=True)
            CurveSketching.plot?
        """
        plot_f = plot(self.f, xmin, xmax, detect_poles='show', ymin=ymin, ymax=ymax, legend_label='f = {}'.format(self.f), legend_color='blue')
        
        for i in self.horizontal_asymptotes():
            if abs(i) != oo:
                plot_f += plot(i, xmin, xmax, ymin=ymin, ymax=ymax, color='lightblue', linestyle='--')
                plot_f += text4plot(' y = {}', i, (xmin, i), col='lightblue')
                
        for i in self.vertical_asymptotes():
            i = round(i, ndigits=2)
            plot_f += text4plot(' x = {}', i, (i, ymin), col = 'lightblue')
        
        for i in self.zeros_solve():
            coord = (i, self.f(x=i))
            plot_f += point(coord, color = 'red', size=10)
            plot_f += text4plot(' {}', round3tuple(coord), coord, col='red', rot = 10)
      
        extrema = self.extrema()
        for i in extrema[0]:
            coord = (i, self.f(x=i))
            plot_f += point(coord, color = 'blue', size=10)
            plot_f += text4plot(' {} Minimum', round3tuple(coord), coord, vert = 'top')
        for i in extrema[1]:
            coord = (i, self.f(x=i))
            plot_f += point(coord, color = 'blue', size=10)
            plot_f += text4plot(' {} Maximum', round3tuple(coord), coord)
        for i in extrema[2]:
            coord = (i, self.f(x=i))
            plot_f += point(coord, color = 'blue', size=10)
            plot_f += text4plot(' {} Inflection Point', round3tuple(coord), coord)
            
        if tan == True:
            tangent = self.tangent()
            plot_f += plot(tangent, x, self.x0 - abs(xmin), self.x0 + xmax, color='green',
                           legend_label='tangent = {}'.format(tangent), legend_color='green') 
        return plot_f
    
    
    def __repr__(self, force_methods = False):
        r"""
        Prints the roots rounded to three digits (solve() method), derivatives, integral, extrema, tangent and asymptotes.
        You can also 'force_methods' i.e. force the roots output from the other three methods (built-in roots(), find_root()
        and self implemented newton).
        
        INPUT:
        
        * "force_methods" - True or False (default: False) - force output of the additional three root solvers
        
        EXAMPLES:
        
            sage: print(CurveSketching((x+17) * (x-3) * (x-1/8)^3).__repr__(force_methods = True))            
            Roots with roots(): [(-17.0, 0.0), (0.125, 0.0), (3.0, 0.0)]
            Roots with find_root(): [(-17.0, 0.0), (3.0, 0.0)]
            Roots with own newton: [(3.0, 0.0)]
            Roots with solve(): [(-17.0, 0.0), (0.125, 0.0), (3.0, 0.0)]
            Derivative: 1/512*(8*x - 1)^3*(x + 17) + 1/512*(8*x - 1)^3*(x - 3) + 3/64*(8*x - 1)^2*(x + 17)*(x - 3)
            2nd Derivative: 1/256*(8*x - 1)^3 + 3/32*(8*x - 1)^2*(x + 17) + 3/32*(8*x - 1)^2*(x - 3) + 3/4*(8*x - 1)*(x + 17)*(x - 3)
            Integral: 1/6*x^6 + 109/40*x^5 - 3597/256*x^4 + 10127/1536*x^3 - 619/512*x^2 + 51/512*x
            Minima: [(2.301, -139.006)]
            Maxima: [(-13.451, 146087.746)]
            Inflection Points: [(0.125, 0.0)]
            Tangent: -619/256*x + 51/512
            Horizontal Limits: [-infinity, +infinity]
            Domain (Vertical Limits): RR\[]
        """
        if self.method == 'solve':
            solve = 'Roots with solve(): {}'.format([round3tuple((i, self.f(x=i))) for i in self.zeros_solve()])
        else:
            solve = 'Roots with roots(): {}'.format([round3tuple((i, self.f(x=i))) for i in self.zeros_roots()])
                
        abl = '\nDerivative: {}'.format(self.derivative())
        abl2 = '\n2nd Derivative: {}'.format(self.derivative(order=2))
        
        inte_check = self.integration()
        if inte_check == None:
            inte = '\nIntegral: {}'.format('Integral is divergent.')
        else:
            inte = '\nIntegral: {}'.format(inte_check)
        
        extrema = self.extrema()
        minmax = '\nMinima: {}\nMaxima: {}\nInflection Points: {}'.format([round3tuple((i, self.f(x=i))) for i in extrema[0]],
                                                                          [round3tuple((i, self.f(x=i))) for i in extrema[1]],
                                                                          [round3tuple((i, self.f(x=i))) for i in extrema[2]])
        tan = '\nTangent: {}'.format(self.tangent())
        hasym = '\nHorizontal Limits: {}'.format([round(i, ndigits=2) for i in self.horizontal_asymptotes()])
        vasym = '\nDomain (Vertical Limits): RR\\{}'.format([round(i, ndigits=2) for i in self.vertical_asymptotes()])
        
        output = solve + abl + abl2 + inte + minmax + tan + hasym + vasym
        
        if force_methods == True:
            try:
                roots = 'Roots with roots(): {}'.format([round3tuple((i, self.f(x=i))) for i in self.zeros_roots()])
            except:
                roots = 'Roots with roots(): {}'.format('Error. Works only on polynomials.')
            try:
                find_root = '\nRoots with find_root(): {}'.format([round3tuple((i, self.f(x=i))) for i in self.zeros_find_root()])
            except:
                find_root = '\nRoots with find_root(): {}'.format('Error')
            try:
                newton = '\nRoots with own newton: {}\n'.format([round3tuple((i, self.f(x=i))) for i in self.zeros_newton()])
            except:
                newton = '\nRoots with own newton: {}\n'.format('Error')
            return roots + find_root + newton + output
        return output