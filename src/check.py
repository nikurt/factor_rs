import sys

sys.setrecursionlimit(10000)

# python3 compatibility
try:
    foo = long
except:
    long = int

ate_loop_count = 29793968203157093288
log_ate_loop_count = 63

curve_order = 21888242871839275222246405745257275088548364400416034343698204186575808495617

# The prime modulus of the field
field_modulus = 21888242871839275222246405745257275088696311157297823662689037894645226208583
# See, it's prime!
assert pow(2, field_modulus, field_modulus) == 2

# The modulus of the polynomial in this representation of FQ12
FQ12_modulus_coeffs = [82, 0, 0, 0, 0, 0, -18, 0, 0, 0, 0, 0]  # Implied + [1]

# Curve order should be prime
assert pow(2, curve_order, curve_order) == 2
# Curve order should be a factor of field_modulus**12 - 1
assert (field_modulus ** 12 - 1) % curve_order == 0


# Elliptic curve doubling
def double(pt):
    x, y = pt
    l1 = (3 * x ** 2)
    l2 = (2 * y)
    l = l1 / l2
    # print("double:: l1 = %d, l2 = %d, l = %d" % (l1.n,l2.n,l.n))
    # py double:: l1 = 7039352064791668362031749933335625902195271755917240662349352063285156878006, l2 = 13110391464550333261117142677974698181990097963069030438958550202690327139681, l = 5255088070853604200184402553820725649001292590477012459392196251894549844818
    # rs double:: l1 = 7039352064791668362031749933335625902195271755917240662349352063285156878006, l2 = 13110391464550333261117142677974698181990097963069030438958550202690327139681, l = 9310243049803878795719365102297529377461416771593574424814731760753039277785

    newx = l ** 2 - 2 * x
    # print("newx",newx)
    newy = -l * newx + l * x - y
    # print("newy",newy)
    # print("double (%d,%d) = (%d,%d)" % (x.n,y.n,newx.n,newy.n))
    return newx, newy


# Elliptic curve addition
def add(p1, p2):
    if p1 is None or p2 is None:
        return p1 if p2 is None else p2
    x1, y1 = p1
    x2, y2 = p2
    if x2 == x1 and y2 == y1:
        return double(p1)
    elif x2 == x1:
        return None
    else:
        l = (y2 - y1) / (x2 - x1)
    newx = l ** 2 - x1 - x2
    newy = -l * newx + l * x1 - y1
    assert newy == (-l * newx + l * x2 - y2)
    return (newx, newy)


# Elliptic curve point multiplication
def multiply(pt, n):
    if n == 0:
        return None
    elif n == 1:
        # print("multiply (%d,%d) * %d = (%d,%d)" % (pt[0].n,pt[1].n,n,pt[0].n,pt[1].n))
        return pt
    elif not n % 2:
        res = multiply(double(pt), n // 2)
        # print("multiply (%d,%d) * %d = (%d,%d)" % (pt[0].n,pt[1].n,n,res[0].n,res[1].n))
        return res
    else:
        res = add(multiply(double(pt), int(n // 2)), pt)
        # print("multiply (%d,%d) * %d = (%d,%d)" % (pt[0].n,pt[1].n,n,res[0].n,res[1].n))
        return res


def eq(p1, p2):
    return p1 == p2


# Extended euclidean algorithm to find modular inverses for
# integers
def inv(a, n):
    if a == 0:
        return 0
    lm, hm = 1, 0
    low, high = a % n, n
    while low > 1:
        r = high // low
        nm, new = hm - lm * r, high - low * r
        lm, low, hm, high = nm, new, lm, low
    res = lm % n
    return res


# A class for field elements in FQ. Wrap a number in this class,
# and it becomes a field element.
class FQ():
    def __init__(self, n):
        if isinstance(n, self.__class__):
            self.n = n.n
        else:
            self.n = n % field_modulus
        assert isinstance(self.n, (int, long))

    def __add__(self, other):
        on = other.n if isinstance(other, FQ) else other
        return FQ((self.n + on) % field_modulus)

    def __mul__(self, other):
        on = other.n if isinstance(other, FQ) else other
        return FQ((self.n * on) % field_modulus)

    def __rmul__(self, other):
        return self * other

    def __radd__(self, other):
        return self + other

    def __rsub__(self, other):
        on = other.n if isinstance(other, FQ) else other
        return FQ((on - self.n) % field_modulus)

    def __sub__(self, other):
        on = other.n if isinstance(other, FQ) else other
        return FQ((self.n - on) % field_modulus)

    def __div__(self, other):
        on = other.n if isinstance(other, FQ) else other
        assert isinstance(on, (int, long))
        r = inv(on, field_modulus)
        res = self.n * r % field_modulus
        # div:: n = 1, r = 2432026985759919469138489527250808343188479017477535962521004210516136245398, res = 2432026985759919469138489527250808343188479017477535962521004210516136245398
        # div:: n = 1, r = 21087453498479301738505683583845423561061080261299122796980902361914303298513, res = 21087453498479301738505683583845423561061080261299122796980902361914303298513
        # div:: n = 2432026985759919469138489527250808343188479017477535962521004210516136245398, r = 21087453498479301738505683583845423561061080261299122796980902361914303298513, res = 14681138511599513868579906292550611339979233093309515871315818100066920017952
        # div:: n = 7039352064791668362031749933335625902195271755917240662349352063285156878006, r = 16560050664832848086065958496307532821520228819545981244151064709290132795008, res = 5255088070853604200184402553820725649001292590477012459392196251894549844818

        # div:: n = Fq(7039352064791668362031749933335625902195271755917240662349352063285156878006), on = 13110391464550333261117142677974698181990097963069030438958550202690327139681, r = 1054657804411666725257106322582614636547841194570646051982474840575437076969, res = Fq(9310243049803878795719365102297529377461416771593574424814731760753039277785)

        # div:: n = Fq(7039352064791668362031749933335625902195271755917240662349352063285156878006), r = 1054657804411666725257106322582614636547841194570646051982474840575437076969, res = Fq(9310243049803878795719365102297529377461416771593574424814731760753039277785)
        # div:: n = Fq(18157769254635907848672757363443528453784178272453828860625883593881387935754), r = 8692195947462038270847387329185208184755276233683606343026879072045979384467, res = Fq(19487853646855803502717319911482627902895396954527750022296372633161786737343)
        # div:: n = Fq(18570508671094861308970420412272049205591995783112888978872489228139002722452), r = 8599263951999741498390012823286631906242527625458707136491257400859549630443, res = Fq(6994108978470069304800728129607443492530055777642244525748264677726031455285)
        # div:: n = Fq(19760575211254260275704618158517966207671483223268603808029365829517623755924), r = 7596631854654444961898591000520960977147355600744279130716907949596604157035, res = Fq(965801065418469171263133164086113828032126717469745510805124821426858127823)
        # div:: n = Fq(14070212845962573227352134915017282957464754993847008570473035153566817473167), r = 2403402005126410256992458165464410967293657589130938487627104004055904400878, res = Fq(16457213427721897696214270963970195665661812259497565867914244341560001147881)
        # div:: n = Fq(14018277940397341683014465763091036605547806946341036878573379237255260306037), r = 2010887891689944622853050607306397980047418624741976773252734516863639070570, res = Fq(968528917167320758575981073797768279419183258954161452223302712869385104095)
        # div:: n = Fq(7167214457385733104678550657788645938365847644176669044945159123786037616001), r = 9188553152785312259936199299562585621724760807205193208373393676046947060214, res = Fq(13703304369351079742367002533642360681342927765371656528473105053759294426850)

        # print("div:: n = %d, on = %d, r = %d, res = %d" % (self.n, on, r, res))
        return FQ(res)

    def __truediv__(self, other):
        return self.__div__(other)

    def __rdiv__(self, other):
        on = other.n if isinstance(other, FQ) else other
        assert isinstance(on, (int, long)), on
        return FQ(inv(self.n, field_modulus) * on % field_modulus)

    def __rtruediv__(self, other):
        return self.__rdiv__(other)

    def __pow__(self, other):
        if other == 0:
            return FQ(1)
        elif other == 1:
            return FQ(self.n)
        elif other % 2 == 0:
            return (self * self) ** (other // 2)
        else:
            return ((self * self) ** int(other // 2)) * self

    def __eq__(self, other):
        if isinstance(other, FQ):
            return self.n == other.n
        else:
            return self.n == other

    def __ne__(self, other):
        return not self == other

    def __neg__(self):
        return FQ(-self.n)

    def __repr__(self):
        return repr(self.n)

    @classmethod
    def one(cls):
        return cls(1)

    @classmethod
    def zero(cls):
        return cls(0)


# Utility methods for polynomial math
def deg(p):
    d = len(p) - 1
    while p[d] == 0 and d:
        d -= 1
    return d


def poly_rounded_div(a, b):
    dega = deg(a)
    degb = deg(b)
    temp = [x for x in a]
    o = [0 for x in a]
    for i in range(dega - degb, -1, -1):
        o[i] += temp[degb + i] / b[degb]
        for c in range(degb + 1):
            temp[c + i] -= o[c]
    return o[:deg(o) + 1]


# A class for elements in polynomial extension fields
class FQP():
    def __init__(self, coeffs, modulus_coeffs):
        assert len(coeffs) == len(modulus_coeffs)
        self.coeffs = [FQ(c) for c in coeffs]
        # The coefficients of the modulus, without the leading [1]
        self.modulus_coeffs = modulus_coeffs
        # The degree of the extension field
        self.degree = len(self.modulus_coeffs)

    def __add__(self, other):
        assert isinstance(other, self.__class__)
        return self.__class__([x + y for x, y in zip(self.coeffs, other.coeffs)])

    def __sub__(self, other):
        assert isinstance(other, self.__class__)
        return self.__class__([x - y for x, y in zip(self.coeffs, other.coeffs)])

    def __mul__(self, other):
        # print("__mul__!!: %s * %s" % (str(self),str(other)))
        if isinstance(other, (FQ, int, long)):
            # print("__mul__!!1")
            return self.__class__([c * other for c in self.coeffs])
        else:
            # print("__mul__!!2")
            assert isinstance(other, self.__class__)
            b = [FQ(0) for i in range(self.degree * 2 - 1)]
            for i in range(self.degree):
                for j in range(self.degree):
                    v = self.coeffs[i] * other.coeffs[j]
                    ind = i + j
                    # print("__mul__!!2, ind: %d v: %s, b[ind]: %s, self.coeffs[i]: %s, other.coeffs[j]: %s" % (ind, str(v), b[ind], self.coeffs[i], other.coeffs[j]))
                    b[ind] += v
                    # b[i + j] += self.coeffs[i] * other.coeffs[j]
            # print("__mul__!!2 b:%s" % b)
            while len(b) > self.degree:
                exp, top = len(b) - self.degree - 1, b.pop()
                for i in range(self.degree):
                    b[exp + i] -= top * FQ(self.modulus_coeffs[i])
            # print("__mul__!!2 b2:%s" % b)
            return self.__class__(b)

    def __rmul__(self, other):
        return self * other

    # Extended euclidean algorithm used to find the modular inverse
    def inv(self):
        lm, hm = [1] + [0] * self.degree, [0] * (self.degree + 1)
        low, high = self.coeffs + [0], self.modulus_coeffs + [1]
        while deg(low):
            r = poly_rounded_div(high, low)
            r += [0] * (self.degree + 1 - len(r))
            nm = [x for x in hm]
            new = [x for x in high]
            assert len(lm) == len(hm) == len(low) == len(high) == len(nm) == len(new) == self.degree + 1
            for i in range(self.degree + 1):
                for j in range(self.degree + 1 - i):
                    nm[i + j] -= lm[i] * r[j]
                    new[i + j] -= low[i] * r[j]
            # print("inv!! lm: %s, low: %s, hm: %s, high: %s, nm: %s, new: %s" % (lm,low,hm,high,nm,new))
            lm, low, hm, high = nm, new, lm, low
        # print("inv!! self: %s low[0]: %s" % (self,low[0]));
        return self.__class__(lm[:self.degree]) / low[0]

    def __div__(self, other):
        if isinstance(other, (FQ, int, long)):
            res = [c / other for c in self.coeffs]
            # print("__div__!!1: self: %s, other: %s, res: %s" % (self,other,res))
            return self.__class__(res)
        else:
            assert isinstance(other, self.__class__)
            return self * other.inv()

    def __truediv__(self, other):
        return self.__div__(other)

    def __pow__(self, other):
        if other == 0:
            # print("pow!!0 %s^%d" % (self,other))
            return self.__class__([1] + [0] * (self.degree - 1))
        elif other == 1:
            # print("pow!!1 %s^%d" % (self,other))
            return self.__class__(self.coeffs)
        elif other % 2 == 0:
            # print("pow!!2 %s^%d" % (self,other))
            return (self * self) ** (other // 2)
        else:
            # print("pow!!3 %s^%d" % (self,other))
            return ((self * self) ** int(other // 2)) * self

    def __repr__(self):
        return repr(self.coeffs)

    def __eq__(self, other):
        assert isinstance(other, self.__class__)
        for c1, c2 in zip(self.coeffs, other.coeffs):
            if c1 != c2:
                return False
        return True

    def __ne__(self, other):
        return not self == other

    def __neg__(self):
        return self.__class__([-c for c in self.coeffs])

    @classmethod
    def one(cls):
        return cls([1] + [0] * (cls.degree - 1))

    @classmethod
    def zero(cls):
        return cls([0] * cls.degree)


# The quadratic extension field
class FQ2(FQP):
    def __init__(self, coeffs):
        self.coeffs = [FQ(c) for c in coeffs]
        self.modulus_coeffs = [1, 0]
        self.degree = 2
        self.__class__.degree = 2


# The 12th-degree extension field
class FQ12(FQP):
    def __init__(self, coeffs):
        self.coeffs = [FQ(c) for c in coeffs]
        self.modulus_coeffs = FQ12_modulus_coeffs
        self.degree = 12
        self.__class__.degree = 12


# Check if a point is the point at infinity
def is_inf(pt):
    return pt is None


# Check that a point is on the curve defined by y**2 == x**3 + b
def is_on_curve(pt, b):
    print("is_in_curve %s %s" % (str(pt), str(b)))
    if is_inf(pt):
        return True
    x, y = pt
    left = y ** 2 - x ** 3
    right = b
    print("is_in_curve %s %s" % (str(left), str(right)))
    return left == right


# Curve is y**2 = x**3 + 3
b = FQ(3)
# Twisted curve over FQ**2
b2 = FQ2([3, 0]) / FQ2([9, 1])
# Extension curve over FQ**12; same b value as over FQ
b12 = FQ12([3] + [0] * 11)

# Generator for curve over FQ
G1 = (FQ(1), FQ(2))
# Generator for twisted curve over FQ2
G2 = (FQ2([10857046999023057135944570762232829481370756359578518086990519993285655852781,
           11559732032986387107991004021392285783925812861821192530917403151452391805634]),
      FQ2([8495653923123431417604973247489272438418190587263600148770280649306958101930,
           4082367875863433681332203403145435568316851327593401208105741076214120093531]))


# Create a function representing the line between P1 and P2,
# and evaluate it at T
def linefunc(P1, P2, T):
    assert P1 and P2 and T  # No points-at-infinity allowed, sorry
    x1, y1 = P1
    x2, y2 = P2
    xt, yt = T
    res = None
    if x1 != x2:
        m = (y2 - y1) / (x2 - x1)
        res = m * (xt - x1) - (yt - y1)
    elif y1 == y2:
        m = 3 * x1 ** 2 / (2 * y1)
        res = m * (xt - x1) - (yt - y1)
    else:
        res = xt - x1
    print("linefunc p1: %s, p2: %s, t: %s, res: %s" % (P1, P2, T, res))
    return res


def cast_point_to_fq12(pt):
    if pt is None:
        return None
    x, y = pt
    return (FQ12([x.n] + [0] * 11), FQ12([y.n] + [0] * 11))


# Check consistency of the "line function"
one, two, three = G1, double(G1), multiply(G1, 3)
negone, negtwo, negthree = multiply(G1, curve_order - 1), multiply(G1, curve_order - 2), multiply(G1, curve_order - 3)

assert linefunc(one, two, one) == FQ(0)
assert linefunc(one, two, two) == FQ(0)
assert linefunc(one, two, three) != FQ(0)
assert linefunc(one, two, negthree) == FQ(0)
assert linefunc(one, negone, one) == FQ(0)
assert linefunc(one, negone, negone) == FQ(0)
assert linefunc(one, negone, two) != FQ(0)
assert linefunc(one, one, one) == FQ(0)
assert linefunc(one, one, two) != FQ(0)
assert linefunc(one, one, negtwo) == FQ(0)


# Main miller loop
def miller_loop(Q, P):
    print("miller_loop: Q: %s, P: %s" % (Q, P))
    if Q is None or P is None:
        return FQ12.one()
    R = Q
    f = FQ12.one()
    print("miller_loop: f: %s" % f)
    for i in range(log_ate_loop_count, -1, -1):
        f = f * f * linefunc(R, R, P)
        R = double(R)
        print("miller_loop: i: %d f: %s, r: %s" % (i, f, R))
        if ate_loop_count & (2 ** i):
            f = f * linefunc(R, Q, P)
            R = add(R, Q)
            print("miller_loop: @i: %s f: %s, r: %s" % (i, f, R))
    # assert R == multiply(Q, ate_loop_count)
    Q1 = (Q[0] ** field_modulus, Q[1] ** field_modulus)
    # assert is_on_curve(Q1, b12)
    nQ2 = (Q1[0] ** field_modulus, -Q1[1] ** field_modulus)
    print("miller_loop q2: %s, q1: %s, nq2: %s, p2: %s, r: %s" % (Q, Q1, nQ2, P, R))
    # assert is_on_curve(nQ2, b12)
    we = linefunc(R, Q1, P)
    f = f * we
    R = add(R, Q1)
    print("miller_loop: we: %s, f: %s, r: %s" % (we, f, R))
    f = f * linefunc(R, nQ2, P)
    print("miller_loop: f: %s, r: %s" % (f, R))
    # R = add(R, nQ2) This line is in many specifications but it technically does nothing
    res = final_exponentiate(f)
    print("miller_loop: res: %s" % res)
    return res


# Pairing computation
def pairing(Q, P):
    print("pairing::\nQ: %s\nP: %s\nb2: %s\nb: %s" % (str(Q), str(P), str(b2), str(b)))
    print("!1")
    assert is_on_curve(P, b)
    print("!2")
    # assert is_on_curve(Q, b2)
    print("!3")
    return miller_loop(twist(Q), cast_point_to_fq12(P))


def final_exponentiate(x):
    p = ((field_modulus ** 12 - 1) // curve_order)
    res = x ** p
    print("final_exponentiate:: f: %s, p: %d, res: %s" % (x, p, res))
    return res


assert is_on_curve(G1, b)
assert is_on_curve(G2, b2)

# "Twist" a point in E(FQ2) into a point in E(FQ12)
w = FQ12([0, 1] + [0] * 10)


# Convert P => -P
def neg(pt):
    if pt is None:
        return None
    x, y = pt
    return (x, -y)


def twist(pt):
    if pt is None:
        return None
    _x, _y = pt
    # Field isomorphism from Z[p] / x**2 to Z[p] / x**2 - 18*x + 82
    xcoeffs = [_x.coeffs[0] - _x.coeffs[1] * 9, _x.coeffs[1]]
    ycoeffs = [_y.coeffs[0] - _y.coeffs[1] * 9, _y.coeffs[1]]
    # Isomorphism into subfield of Z[p] / w**12 - 18 * w**6 + 82,
    # where w**6 = x
    nx = FQ12([xcoeffs[0]] + [0] * 5 + [xcoeffs[1]] + [0] * 5)
    ny = FQ12([ycoeffs[0]] + [0] * 5 + [ycoeffs[1]] + [0] * 5)
    # Divide x coord by w**2 and y coord by w**3
    return (nx * w ** 2, ny * w ** 3)


G12 = twist(G2)
# Check that the twist creates a point that is on the curve
assert is_on_curve(G12, b12)

"""
# A class for elements in polynomial extension fields
class FQP():
    def __init__(self, coeffs, modulus_coeffs):
        assert len(coeffs) == len(modulus_coeffs)
        self.coeffs = [FQ(c) for c in coeffs]
        # The coefficients of the modulus, without the leading [1]
        self.modulus_coeffs = modulus_coeffs
        # The degree of the extension field
        self.degree = len(self.modulus_coeffs)

    def __add__(self, other):
        assert isinstance(other, self.__class__)
        return self.__class__([x+y for x,y in zip(self.coeffs, other.coeffs)])

    def __sub__(self, other):
        assert isinstance(other, self.__class__)
        return self.__class__([x-y for x,y in zip(self.coeffs, other.coeffs)])

    def __mul__(self, other):
        if isinstance(other, (FQ, int, long)):
            return self.__class__([c * other for c in self.coeffs])
        else:
            assert isinstance(other, self.__class__)
            b = [FQ(0) for i in range(self.degree * 2 - 1)]
            for i in range(self.degree):
                for j in range(self.degree):
                    b[i + j] += self.coeffs[i] * other.coeffs[j]
            while len(b) > self.degree:
                exp, top = len(b) - self.degree - 1, b.pop()
                for i in range(self.degree):
                    b[exp + i] -= top * FQ(self.modulus_coeffs[i])
            return self.__class__(b)

    def __rmul__(self, other):
        return self * other

    def __div__(self, other):
        if isinstance(other, (FQ, int, long)):
            return self.__class__([c / other for c in self.coeffs])
        else:
            assert isinstance(other, self.__class__)
            return self * other.inv()

    def __truediv__(self, other):
        return self.__div__(other)

    def __pow__(self, other):
        if other == 0:
            return self.__class__([1] + [0] * (self.degree - 1))
        elif other == 1:
            return self.__class__(self.coeffs)
        elif other % 2 == 0:
            return (self * self) ** (other // 2)
        else:
            return ((self * self) ** int(other // 2)) * self

    # Extended euclidean algorithm used to find the modular inverse
    def inv(self):
        lm, hm = [1] + [0] * self.degree, [0] * (self.degree + 1)
        low, high = self.coeffs + [0], self.modulus_coeffs + [1]
        while deg(low):
            r = poly_rounded_div(high, low)
            r += [0] * (self.degree + 1 - len(r))
            nm = [x for x in hm]
            new = [x for x in high]
            assert len(lm) == len(hm) == len(low) == len(high) == len(nm) == len(new) == self.degree + 1
            for i in range(self.degree + 1):
                for j in range(self.degree + 1 - i):
                    nm[i+j] -= lm[i] * r[j]
                    new[i+j] -= low[i] * r[j]
            lm, low, hm, high = nm, new, lm, low
        return self.__class__(lm[:self.degree]) / low[0]

    def __repr__(self):
        return repr(self.coeffs)

    def __eq__(self, other):
        assert isinstance(other, self.__class__)
        for c1, c2 in zip(self.coeffs, other.coeffs):
            if c1 != c2:
                return False
        return True

    def __ne__(self, other):
        return not self == other

    def __neg__(self):
        return self.__class__([-c for c in self.coeffs])

    @classmethod
    def one(cls):
        return cls([1] + [0] * (cls.degree - 1))

    @classmethod
    def zero(cls):
        return cls([0] * cls.degree)

# The quadratic extension field
class FQ2(FQP):
    def __init__(self, coeffs):
        self.coeffs = [FQ(c) for c in coeffs]
        self.modulus_coeffs = [1, 0]
        self.degree = 2
        self.__class__.degree = 2

# The 12th-degree extension field
class FQ12(FQP):
    def __init__(self, coeffs):
        self.coeffs = [FQ(c) for c in coeffs]
        self.modulus_coeffs = FQ12_modulus_coeffs
        self.degree = 12
        self.__class__.degree = 12
"""

# p1:[G1Point { x: Fq(14381381155403816334798493340013261768598472201420748168820652964797987547462), y: Fq(18869406103595735495698386560184662301160822785256257928525476998858093794620) }, G1Point { x: Fq(12069817318037655725670852344933769708171796463120307972343537074143694270596), y: Fq(3231429222134733375287079359249040988942692960988404421142202529154679566108) }, G1Point { x: Fq(21775469447366302559279280857827659314936507430950512123672221590066690909025), y: Fq(17656487425289707757904355540441031854480107147912420173634687575335079425665) }, G1Point { x: Fq(12411515119042684284347439080261290708250768953338792298174140687146844441182), y: Fq(20138364572167484814256846272295865295005456643842414091875980986217735719454) }]
# p2:[G2Point { x: [Fq(1222372640101605630716805870247697620367955331142053227617023362218224623675), Fq(20145010989865280903946427053568635653229019267312571214695898029471466643189)], y: [Fq(18627856813852125933293899443711247288107709517554840445817048145846492416046), Fq(7179793276317950049211890536373035911436265064183990876740337674144380042597)] }, G2Point { x: [Fq(6138690762340449119094879065809607641134675188862351822918120352834699055236), Fq(19891240116629936231754768188959733167757134087113133540973038412946691553537)], y: [Fq(13005312490757206815920942747588625201931220927023222301565124660940116540363), Fq(9652081022470897254611722611899202761684115862349335883608737443091937419563)] }, G2Point { x: [Fq(11559732032986387107991004021392285783925812861821192530917403151452391805634), Fq(10857046999023057135944570762232829481370756359578518086990519993285655852781)], y: [Fq(4082367875863433681332203403145435568316851327593401208105741076214120093531), Fq(8495653923123431417604973247489272438418190587263600148770280649306958101930)] }, G2Point { x: [Fq(14300072486619654782841429113650020748800347468876509540346977848520405442186), Fq(13034590777758516932550084806338576045647608732831270584871813280920866673420)], y: [Fq(5846437282942158081535633083579968858933740559240158365995208260434325310477), Fq(4326610515916630131482228944454373446417733591398164445981042744194140563287)] }]

p1 = [
    (FQ(14381381155403816334798493340013261768598472201420748168820652964797987547462),
     FQ(18869406103595735495698386560184662301160822785256257928525476998858093794620)),
    (FQ(12069817318037655725670852344933769708171796463120307972343537074143694270596),
     FQ(3231429222134733375287079359249040988942692960988404421142202529154679566108)),
    (FQ(4175458037543491544904263421106282279888787737666823543403682962452691369986),
     FQ(10396040821576636452323344393044999842161477672879106129746340283771235116737)),
    (FQ(12411515119042684284347439080261290708250768953338792298174140687146844441182),
     FQ(20138364572167484814256846272295865295005456643842414091875980986217735719454)),
]
p2 = [
    (FQ2([1222372640101605630716805870247697620367955331142053227617023362218224623675,
          20145010989865280903946427053568635653229019267312571214695898029471466643189]),
     FQ2([18627856813852125933293899443711247288107709517554840445817048145846492416046,
          7179793276317950049211890536373035911436265064183990876740337674144380042597])),
    (FQ2([6138690762340449119094879065809607641134675188862351822918120352834699055236,
          19891240116629936231754768188959733167757134087113133540973038412946691553537]),
     FQ2([13005312490757206815920942747588625201931220927023222301565124660940116540363,
          9652081022470897254611722611899202761684115862349335883608737443091937419563])),
    (FQ2([11559732032986387107991004021392285783925812861821192530917403151452391805634,
          10857046999023057135944570762232829481370756359578518086990519993285655852781]),
     FQ2([4082367875863433681332203403145435568316851327593401208105741076214120093531,
          8495653923123431417604973247489272438418190587263600148770280649306958101930])),
    (FQ2([14300072486619654782841429113650020748800347468876509540346977848520405442186,
          13034590777758516932550084806338576045647608732831270584871813280920866673420]),
     FQ2([5846437282942158081535633083579968858933740559240158365995208260434325310477,
          4326610515916630131482228944454373446417733591398164445981042744194140563287])),
]

print("----------------")


# def ff():
#    Q=(FQ12( [FQ(0), FQ(0), FQ(6429969845989378501122694582889211781453424659884192589945847793458736424321), FQ(0), FQ(0), FQ(0), FQ(0), FQ(0), FQ(13034590777758516932550084806338576045647608732831270584871813280920866673420), FQ(0), FQ(0), FQ(0)] ), FQ12( [FQ(0), FQ(0), FQ(0), FQ(10683428383371037342688384074005158018566760551252325677543899351977512658060), FQ(0), FQ(0), FQ(0), FQ(0), FQ(0), FQ(4326610515916630131482228944454373446417733591398164445981042744194140563287), FQ(0), FQ(0)] ))
#    P=(FQ12([FQ(12411515119042684284347439080261290708250768953338792298174140687146844441182), FQ(0), FQ(0), FQ(0), FQ(0), FQ(0), FQ(0), FQ(0), FQ(0), FQ(0), FQ(0), FQ(0)] ), FQ12([FQ(20138364572167484814256846272295865295005456643842414091875980986217735719454), FQ(0), FQ(0), FQ(0), FQ(0), FQ(0), FQ(0), FQ(0), FQ(0), FQ(0), FQ(0), FQ(0)]))
#    res = miller_loop(Q,P)
#    print(res)

# ff()
# print("----------------")
# inv!! lm: [Fq(1), Fq(0), Fq(0)], low: [Fq(9), Fq(1), Fq(0)], hm: [Fq(0), Fq(0), Fq(0)], high: [Fq(1), Fq(0), Fq(1)], nm: [Fq(0), Fq(21888242871839275222246405745257275088696311157297823662689037894645226208582), Fq(0)], new: [Fq(1), Fq(21888242871839275222246405745257275088696311157297823662689037894645226208574), Fq(0)]
# inv!! lm: [Fq(0), Fq(21888242871839275222246405745257275088696311157297823662689037894645226208582), Fq(0)], low: [Fq(1), Fq(21888242871839275222246405745257275088696311157297823662689037894645226208574), Fq(0)], hm: [Fq(1), Fq(0), Fq(0)], high: [Fq(9), Fq(1), Fq(0)], nm: [Fq(1), Fq(2432026985759919469138489527250808343188479017477535962521004210516136245398), Fq(0)], new: [Fq(19456215886079355753107916218006466745507832139820287700168033684129089963194), Fq(0), Fq(0)]
# inv!! self: FQP { coeffs: [Fqpairing(9), Fq(1)], modulus_coeffs: [Fq(1), Fq(0)], degree: 2 } low[0]: Fq(19456215886079355753107916218006466745507832139820287700168033684129089963194)


# print(b2)

def doPairing(p1, p2):
    print("doPairing:\np1:%s\np2:%s" % (str(p1), str(p2)))
    assert len(p1) == len(p2)
    product = FQ12.one()
    for i in range(len(p1)):
        zr = pairing(p2[i], p1[i])
        rz = product * zr
        print("pairing:: product_%d: %s, zr: %s, rz: %s" % (i, product, zr, rz))
        product = rz
    return product.coeffs[0] == 1


doPairing(p1, p2)
