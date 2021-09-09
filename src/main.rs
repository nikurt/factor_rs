use std::convert::TryFrom;
use std::fmt;
use core::ops::Shr;
use primitive_types::U256;
use primitive_types::U512;
use std::ops::Div;
use num_bigint::BigUint;
use num_traits::{Num,One,Zero};
use num_traits::ops::checked::CheckedSub;
use num_integer::Integer;

// Curve is y**2 = x**3 + 3
fn b() -> Fq {
    Fq::new(U256::from(3))
}

// Twisted curve over FQ**2
fn b2() -> FQP {
    FQP::fq2(&[3, 0]).div_fqp(FQP::fq2(&[9, 1]))
}

// Extension curve over FQ**12; same b value as over FQ
fn b12() -> FQP {
    FQP::fq12(&[3, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0])
}

#[derive(Eq, PartialEq, Debug, Default, Clone, Copy)]
struct Fq(U256);

// TODO: Use BigInt instead.
#[derive(Copy, Clone, Eq, PartialEq, Debug)]
struct U256WithSign {
    x: U256,
    is_neg: bool,
}

impl fmt::Display for U256WithSign {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        if self.is_neg {
            write!(f, "-{}", self.x)
        } else {
            write!(f, "{}", self.x)
        }
    }
}

impl U256WithSign {
    fn sub(self, x: U256WithSign) -> U256WithSign {
        if x.x.is_zero() {
            self
        } else if self == x {
            U256WithSign { x: U256::zero(), is_neg: false }
        } else if self.is_neg == x.is_neg {
            if self.x >= x.x {
                U256WithSign { x: self.x - x.x, is_neg: x.is_neg }
            } else {
                U256WithSign { x: self.x.checked_add(x.x).unwrap(), is_neg: !self.is_neg }
            }
        } else {
            U256WithSign { x: self.x.checked_add(x.x).unwrap(), is_neg: self.is_neg }
        }
    }

    fn mul(self, x: U256WithSign) -> U256WithSign {
        if self.x.is_zero() || x.x.is_zero() {
            U256WithSign { x: U256::zero(), is_neg: false }
        } else {
            U256WithSign { x: self.x * x.x, is_neg: self.is_neg ^ x.is_neg }
        }
    }

    fn div(self, x: U256WithSign) -> U256WithSign {
        let (a, b) = self.x.div_mod(x.x);
        if b.is_zero() {
            if a.is_zero() {
                U256WithSign { x: U256::zero(), is_neg: false }
            } else {
                U256WithSign { x: a, is_neg: self.is_neg ^ x.is_neg }
            }
        } else {
            let is_neg = self.is_neg ^ x.is_neg;
            if !is_neg {
                U256WithSign { x: a, is_neg: false }
            } else {
                U256WithSign { x: a + U256::one(), is_neg: is_neg }
            }
        }
    }
}


impl fmt::Display for Fq {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        write!(f, "Fq({})", self.0)
    }
}

fn inv(a: U256, b: U256) -> U256 {
    if a.is_zero() {
        return a;
    }
    /*
    let mut x = U256::one();
    let mut y = U256::zero();
    let mut x1 = U256::zero();
    let mut y1 = U256::one();
    let mut a1 = a;
    let mut b1 = b;

    while b1 > U256::zero() {
        let q = a1.div(b1);
        println!("x: {}, q: {}, x1: {}",x,q,x1);
        let nx = x-q*x1;
        x = x1;
        x1=nx;
        println!("y: {}, q: {}, y1: {}",y,q,y1);
        let ny=y-q*y1;
        y = y1;
        y1 = ny;
        println!("a1: {}, q: {}, b1: {}",a1,q,b1);
        let nab = a1-q*b1;
        a1=b1;
        b1=nab;
    }
    a1
     */


    let mut lm = U256WithSign { x: U256::one(), is_neg: false };
    let mut hm = U256WithSign { x: U256::zero(), is_neg: false };
    let mut low = U256WithSign { x: a % b, is_neg: false };
    let mut high = U256WithSign { x: b, is_neg: false };
    println!("inv:: a: {}, n: {}, lm: {}, hm: {}, low: {}, high: {}", a, b, lm, hm, low, high);
    loop {
        let c = low.sub(U256WithSign { x: U256::one(), is_neg: false });
        if c.is_neg || c.x.is_zero() {
            break;
        }
        let r = high.div(low);
        // println!("lm: {:?}, r: {:?}, hm: {:?}",lm,r,hm);
        let nm = hm.sub(lm.mul(r));
        let new = high.sub(low.mul(r));

        high = low;
        low = new;
        hm = lm;
        lm = nm;
        println!("inv:: lm: {}, hm: {}, low: {}, high: {}", lm, hm, low, high);
    }
    let res =
        if lm.is_neg {
            (b - lm.x.div_mod(b).1).div_mod(b).1
        } else {
            lm.x.div_mod(b).1
        };
    println!("inv:: a: {}, n: {}, res: {}", a, b, res);
    res
}

impl Fq {
    fn FIELD_MODULUS() -> U256 {
        let X: U256 = U256::from_dec_str("21888242871839275222246405745257275088696311157297823662689037894645226208583").unwrap();
        X
    }

    fn new(x: U256) -> Self {
        Self(x.div_mod(Fq::FIELD_MODULUS()).1)
    }
    fn new_512(x: U512) -> Self {
        Self(U256::try_from(x.div_mod(U512::from(Fq::FIELD_MODULUS())).1).unwrap())
    }
    fn add(self, x: Fq) -> Fq {
        Fq::new(self.0 + x.0)
    }
    fn sub(self, x: Fq) -> Fq {
        Fq::new(self.0.checked_add(x.neg().0).unwrap())
    }
    fn mul(self, x: Fq) -> Fq {
        Fq::new_512(self.0.full_mul(x.0))
    }
    fn neg(self) -> Fq {
        Fq::new(Fq::FIELD_MODULUS() - self.0)
    }
    fn inv(self) -> Fq {
        Fq::new(inv(self.0, Fq::FIELD_MODULUS()))
    }
    fn div(self, x: Fq) -> Fq {
        let r = inv(x.0, Fq::FIELD_MODULUS());
        let res = Fq::new_512(self.0.full_mul(inv(x.0, Fq::FIELD_MODULUS())));
        println!("div:: n = {}, on = {}, r = {}, res = {}", self, x.0, r, res);
        res
    }
    fn is_zero(&self) -> bool {
        self.0.is_zero()
    }
}

#[derive(Copy, Clone, Eq, PartialEq, Debug)]
struct G1Point {
    x: Fq,
    y: Fq,
}

impl fmt::Display for G1Point {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        write!(f, "G1Point({},{})", self.x, self.y)
    }
}

impl G1Point {
    fn double(self) -> Self {
        println!("double: {}", self);
        let l1 = self.x.mul(self.x).mul(Fq::new(U256::from(3)));
        let l2 = self.y.mul(Fq::new(U256::from(2)));
        let l = l1.div(l2);
        println!("double:: l1 = {}, l2 = {}, l = {}", l1, l2, l);
        let new_x = l.mul(l).sub(self.x.mul(Fq::new(U256::from(2))));
        let new_y = l.neg().mul(new_x).add(l.mul(self.x)).sub(self.y);
        let res = G1Point { x: new_x, y: new_y };
        println!("double: {} = {}", self, res);
        res
    }

    fn add(self, p2: G1Point) -> G1Point {
        println!("add: {} + {}", self, p2);
        let res = if self.x == p2.x && self.y == p2.y {
            self.double()
        } else {
            assert_ne!(self.x, p2.x);
            let l1 = p2.y.sub(self.y);
            let l2 = p2.x.sub(self.x);
            let l = l1.div(l2);
            println!("l1={} l2={} l={}", l1, l2, l);
            let new_x = l.mul(l).sub(self.x).sub(p2.x);
            let new_y = l.neg().mul(new_x).add(l.mul(self.x)).sub(self.y);
            println!("add:: self: {}, p2: {}, l: {}, new_x: {}, new_y: {}", self, p2, l, new_x, new_y);
            assert_eq!(new_y, l.neg().mul(new_x).add(l.mul(p2.x)).sub(p2.y));
            G1Point { x: new_x, y: new_y }
        };
        println!("add: {} + {} = {}", self, p2, res);
        res
    }

    /// @return r the negation of p, i.e. p.addition(p.negate()) should be zero.
    fn neg(self) -> G1Point {
        println!("neg: {}", self);
        // The prime q in the base field F_q for G1
        let res = if self.x.0.is_zero() && self.y.0.is_zero() {
            G1Point {
                x: Fq::new(U256::from(0)),
                y: Fq::new(U256::from(0)),
            }
        } else {
            G1Point {
                x: self.x,
                y: self.y.neg(),
            }
        };
        println!("neg: {} = {}", self, res);
        res
    }

    /// @return r the product of a point on G1 and a scalar, i.e.
    /// p == p.scalar_mul(1) and p.addition(p) == p.scalar_mul(2) for all points p.
    fn mul_fq(self, n: Fq) -> G1Point {
        println!("mul_fq: {}*{}", self, n);
        assert_ne!(n.0, U256::zero());
        let res = if n.0 == U256::one() {
            return self;
        } else if (n.0 % 2).is_zero() {
            self.double().mul_fq(Fq::new(n.0.shr(1)))
        } else {
            self.double().mul_fq(Fq::new(n.0.shr(1))).add(self)
        };
        println!("mul_fq: {}*{} = {}", self, n, res);
        res
    }
}

#[derive(Copy, Clone, Eq, PartialEq, Debug)]
struct G2Point {
    x: [Fq; 2],
    y: [Fq; 2],
}

impl fmt::Display for G2Point {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        write!(f, "G2Point([{},{}],[{},{}])", self.x[0], self.x[1], self.y[0], self.y[1])
    }
}


/// @return the generator of G1
fn P1() -> G1Point {
    G1Point {
        x: Fq::new(U256::from(1)),
        y: Fq::new(U256::from(2)),
    }
}

/// @return the generator of G2
fn P2() -> G2Point {
    // Original code point
    G2Point {
        x: [
            Fq::new(U256::from_dec_str(
                "10857046999023057135944570762232829481370756359578518086990519993285655852781",
            )
                .unwrap()),
            Fq::new(U256::from_dec_str(
                "11559732032986387107991004021392285783925812861821192530917403151452391805634",
            )
                .unwrap()),
        ],
        y: [
            Fq::new(U256::from_dec_str(
                "8495653923123431417604973247489272438418190587263600148770280649306958101930",
            )
                .unwrap()),
            Fq::new(U256::from_dec_str(
                "4082367875863433681332203403145435568316851327593401208105741076214120093531",
            )
                .unwrap()),
        ],
    }
}

fn is_on_curve_g1(p: G1Point, b: Fq) -> bool {
    return p.y.mul(p.y).sub(p.x.mul(p.x).mul(p.x)) == b;
}

fn is_on_curve_g2(p: G2Point, b: FQP) -> bool {
    assert_eq!(b.degree, 2);
    let px = FQP::fq2_g2(p.x);
    let py = FQP::fq2_g2(p.y);
    return *(&py.mul_fqp(&py).sub(&px.mul_fqp(&px).mul_fqp(&px))) == b;
}

fn w() -> FQP {
    FQP::fq12(&[0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0])
}

fn twist(q: G2Point) -> G12Point {
    /*
    _x, _y = pt
    # Field isomorphism from Z[p] / x**2 to Z[p] / x**2 - 18*x + 82
    xcoeffs = [_x.coeffs[0] - _x.coeffs[1] * 9, _x.coeffs[1]]
    ycoeffs = [_y.coeffs[0] - _y.coeffs[1] * 9, _y.coeffs[1]]
    */
    let xcoeffs = [q.x[0].sub(q.x[1].mul(Fq::new(U256::from(9)))), q.x[1]];
    let ycoeffs = [q.y[0].sub(q.y[1].mul(Fq::new(U256::from(9)))), q.y[1]];

    /*
    # Isomorphism into subfield of Z[p] / w**12 - 18 * w**6 + 82,
    # where w**6 = x
    nx = FQ12([xcoeffs[0]] + [0] * 5 + [xcoeffs[1]] + [0] * 5)
    ny = FQ12([ycoeffs[0]] + [0] * 5 + [ycoeffs[1]] + [0] * 5)
     */
    let nx = FQP::fq12_fq(&[xcoeffs[0], Fq::new(U256::zero()), Fq::new(U256::zero()), Fq::new(U256::zero()), Fq::new(U256::zero()), Fq::new(U256::zero()), xcoeffs[1], Fq::new(U256::zero()), Fq::new(U256::zero()), Fq::new(U256::zero()), Fq::new(U256::zero()), Fq::new(U256::zero())]);
    let ny = FQP::fq12_fq(&[ycoeffs[1], Fq::new(U256::zero()), Fq::new(U256::zero()), Fq::new(U256::zero()), Fq::new(U256::zero()), Fq::new(U256::zero()), ycoeffs[1], Fq::new(U256::zero()), Fq::new(U256::zero()), Fq::new(U256::zero()), Fq::new(U256::zero()), Fq::new(U256::zero())]);
    // Divide x coord by w**2 and y coord by w**3
    let res = G12Point { x: nx.mul_fqp(&w()).mul_fqp(&w()), y: ny.mul_fqp(&w()).mul_fqp(&w()).mul_fqp(&w()) };
    res
}

#[derive(Debug, Clone, Eq, PartialEq)]
struct G12Point {
    x: FQP,
    y: FQP,
}

impl G12Point {
    fn double(self) -> Self {
        let l1 = self.x.mul_fqp(&self.x).mul_fq(Fq::new(U256::from(3)));
        let l2 = self.y.mul_fq(Fq::new(U256::from(2)));
        let l = l1.div_fqp(l2);
        let new_x = l.mul_fqp(&l).sub(&self.x.mul_fq(Fq::new(U256::from(2))));
        let new_y = l.mul_fqp(&self.x).sub(&self.y).sub(&l.mul_fqp(&new_x));
        let res = G12Point { x: new_x, y: new_y };
        res
    }

    fn add(self, p2: &G12Point) -> G12Point {
        let res = if self.x == p2.x && self.y == p2.y {
            self.double()
        } else {
            let l1 = p2.y.sub(&self.y);
            let l2 = p2.x.sub(&self.x);
            let l = l1.div_fqp(l2);
            let new_x = l.mul_fqp(&l).sub(&self.x).sub(&p2.x);
            let new_y = l.mul_fqp(&self.x).sub(&self.y).sub(&l.mul_fqp(&new_x));
            G12Point { x: new_x, y: new_y }
        };
        res
    }
}

fn cast_point_to_fq12(p: G1Point) -> G12Point {
    G12Point {
        x: FQP::fq12_fq(&[p.x, Fq::new(U256::zero()), Fq::new(U256::zero()), Fq::new(U256::zero()), Fq::new(U256::zero()), Fq::new(U256::zero()), Fq::new(U256::zero()), Fq::new(U256::zero()), Fq::new(U256::zero()), Fq::new(U256::zero()), Fq::new(U256::zero()), Fq::new(U256::zero())]),
        y: FQP::fq12_fq(&[p.y, Fq::new(U256::zero()), Fq::new(U256::zero()), Fq::new(U256::zero()), Fq::new(U256::zero()), Fq::new(U256::zero()), Fq::new(U256::zero()), Fq::new(U256::zero()), Fq::new(U256::zero()), Fq::new(U256::zero()), Fq::new(U256::zero()), Fq::new(U256::zero())]),
    }
}

fn linefunc(p1: &G12Point, p2: &G12Point, t: &G12Point) -> FQP {
    if p1.x != p2.x {
        let m = p2.y.sub(&p1.y).div_fqp(p2.x.sub(&p1.x));
        m.mul_fqp(&t.x.sub(&p1.x)).sub(&t.y.sub(&p1.y))
    } else if p1.y == p2.y {
        let m = p1.x.mul_fqp(&p1.x).mul_fq(Fq::new(U256::from(3))).div_fqp(p1.y.mul_fq(Fq::new(U256::from(2))));
        m.mul_fqp(&t.x.sub(&p1.x)).sub(&t.y.sub(&p1.y))
    } else {
        t.x.sub(&p1.x)
    }
}


fn miller_loop(q2: G12Point, p2: G12Point) -> FQP {
    let mut r:G12Point = q2.clone();
    let mut f = FQP::fq12_one();
    let ate_loop_count = U256::from_dec_str("29793968203157093288").unwrap();
    const log_ate_loop_count: u64 = 63;
    for i in (0..(log_ate_loop_count + 1)).rev() {
        f = f.mul_fqp(&f).mul_fqp(&linefunc(&r, &r, &p2));
        r = r.double();
        if !(ate_loop_count & U256::from(1 << i)).is_zero() {
            f = f.mul_fqp(&linefunc(&r, &q2, &p2));
            r = r.add(&q2);
        }
    }
    let q1 = G12Point { x: q2.x.pow(Fq::FIELD_MODULUS()), y: q2.y.pow(Fq::FIELD_MODULUS()) };
    let nq2 = G12Point { x: q1.x.pow(Fq::FIELD_MODULUS()), y: (q1.y.neg()).pow(Fq::FIELD_MODULUS()) };
    f = f.mul_fqp(&linefunc(&r, &q1, &p2));
    r = r.add(&q1);
    f = f.mul_fqp(&linefunc(&r, &nq2, &p2));
    final_exponentiate(&f)
}

fn curve_order() -> U256 {
    return U256::from_dec_str("21888242871839275222246405745257275088548364400416034343698204186575808495617").unwrap();
}

fn final_exponentiate(f: &FQP) -> FQP {
    return f.pow_biguint(&BigUint::parse_bytes(Fq::FIELD_MODULUS().to_string().as_bytes(), 10).unwrap().pow(12).checked_sub(&BigUint::one()).unwrap());//.(FIELD_ORDER));
}

fn pairing_p(p: G1Point, q: G2Point) -> FQP {
    assert!(is_on_curve_g1(p, b()));
    assert!(is_on_curve_g2(q, b2()));
    let q2 = twist(q);
    let p2 = cast_point_to_fq12(p);
    FQP::fq2(&[0, 0])
    //let res = miller_loop(q2, p2);
    //res
}

/// @return the result of computing the pairing check
/// e(p1[0], p2[0]) *  .... * e(p1[n], p2[n]) == 1
/// For example pairing([P1(), P1().negate()], [P2(), P2()]) should
/// return true.
fn pairing(p1: &[G1Point], p2: &[G2Point]) -> bool {
    assert!(p1.len() == p2.len());
    let mut product = FQP::fq12(&[1, 0]);
    for i in 0..p1.len() {
        let p = pairing_p(p1[i], p2[i]);
        product = product.mul_fqp(&p);
    }
    return product.coeffs[0].0 == U256::one();
}

/// Convenience method for a pairing check for four pairs.
fn pairingProd4(
    a1: G1Point,
    a2: G2Point,
    b1: G1Point,
    b2: G2Point,
    c1: G1Point,
    c2: G2Point,
    d1: G1Point,
    d2: G2Point,
) -> bool {
    pairing(&[a1, b1, c1, d1], &[a2, b2, c2, d2])
}

struct VerifyingKey {
    alfa1: G1Point,
    beta2: G2Point,
    gamma2: G2Point,
    delta2: G2Point,
    IC: [G1Point; 2],
}

struct Proof {
    A: G1Point,
    B: G2Point,
    C: G1Point,
}

fn verifyingKey() -> VerifyingKey {
    VerifyingKey {
        alfa1: G1Point {
            x: Fq::new(U256::from_dec_str(
                "12069817318037655725670852344933769708171796463120307972343537074143694270596",
            )
                .unwrap()),
            y: Fq::new(U256::from_dec_str(
                "3231429222134733375287079359249040988942692960988404421142202529154679566108",
            )
                .unwrap()),
        },
        beta2: G2Point {
            x: [
                Fq::new(U256::from_dec_str(
                    "6138690762340449119094879065809607641134675188862351822918120352834699055236",
                )
                    .unwrap()),
                Fq::new(U256::from_dec_str(
                    "19891240116629936231754768188959733167757134087113133540973038412946691553537",
                )
                    .unwrap()),
            ],
            y: [
                Fq::new(U256::from_dec_str(
                    "13005312490757206815920942747588625201931220927023222301565124660940116540363",
                )
                    .unwrap()),
                Fq::new(U256::from_dec_str(
                    "9652081022470897254611722611899202761684115862349335883608737443091937419563",
                )
                    .unwrap()),
            ],
        },
        gamma2: G2Point {
            x: [
                Fq::new(U256::from_dec_str(
                    "11559732032986387107991004021392285783925812861821192530917403151452391805634",
                )
                    .unwrap()),
                Fq::new(U256::from_dec_str(
                    "10857046999023057135944570762232829481370756359578518086990519993285655852781",
                )
                    .unwrap()),
            ],
            y: [
                Fq::new(U256::from_dec_str(
                    "4082367875863433681332203403145435568316851327593401208105741076214120093531",
                )
                    .unwrap()),
                Fq::new(U256::from_dec_str(
                    "8495653923123431417604973247489272438418190587263600148770280649306958101930",
                )
                    .unwrap()),
            ],
        },
        delta2: G2Point {
            x: [
                Fq::new(U256::from_dec_str(
                    "14300072486619654782841429113650020748800347468876509540346977848520405442186",
                )
                    .unwrap()),
                Fq::new(U256::from_dec_str(
                    "13034590777758516932550084806338576045647608732831270584871813280920866673420",
                )
                    .unwrap()),
            ],
            y: [
                Fq::new(U256::from_dec_str(
                    "5846437282942158081535633083579968858933740559240158365995208260434325310477",
                )
                    .unwrap()),
                Fq::new(U256::from_dec_str(
                    "4326610515916630131482228944454373446417733591398164445981042744194140563287",
                )
                    .unwrap()),
            ],
        },
        IC: [
            G1Point {
                x: Fq::new(U256::from_dec_str(
                    "11235677443219608181162325740620864972713426694771199166344399448217229053842",
                )
                    .unwrap()),
                y: Fq::new(U256::from_dec_str(
                    "1447324764906200183505263375815319707129765120430582581937303549365705732219",
                )
                    .unwrap()),
            },
            G1Point {
                x: Fq::new(U256::from_dec_str(
                    "18107605096562723854700666435595314810062186832481469946731400411002744132475",
                )
                    .unwrap()),
                y: Fq::new(U256::from_dec_str(
                    "17499317168194804241681774211615986635343204560183427050823794048667776674132",
                )
                    .unwrap()),
            },
        ],
    }
}

fn verify(input: &[U256], proof: &Proof) -> bool {
    let vk = verifyingKey();
    assert!(input.len() + 1 == vk.IC.len());
    // Compute the linear combination vk_x
    let mut vk_x = G1Point { x: Fq::new(U256::from(0)), y: Fq::new(U256::from(0)) };
    println!("vk_x: {},{}", vk_x.x, vk_x.y);
    for i in 0..input.len() {
        assert!(input[i] < Fq::FIELD_MODULUS());
        vk_x = vk_x.add(vk.IC[i + 1].mul_fq(Fq::new(input[i])));
        println!("vk_x: {},{}", vk_x.x, vk_x.y);
    }
    vk_x = vk_x.add(vk.IC[0]);
    println!("vk_x: {},{}", vk_x.x, vk_x.y);
    pairingProd4(
        proof.A.neg(),
        proof.B,
        vk.alfa1,
        vk.beta2,
        vk_x,
        vk.gamma2,
        proof.C,
        vk.delta2,
    )
}

pub struct Verifier {}

impl Verifier {
    /// @return r  bool true if proof is valid
    fn verifyProof(a: [Fq; 2], b: [[Fq; 2]; 2], c: [Fq; 2], input: &[U256]) -> bool {
        let proof = Proof {
            A: G1Point { x: a[0], y: a[1] },
            B: G2Point {
                x: [b[0][0], b[0][1]],
                y: [b[1][0], b[1][1]],
            },
            C: G1Point { x: c[0], y: c[1] },
        };
        verify(input, &proof)
    }
}

fn fq12_modulus_coeffs() -> Vec<U256> {
    vec![82, 0, 0, 0, 0, 0, -18, 0, 0, 0, 0, 0].iter().map(|x| U256::from(*x)).collect()
}

#[derive(Debug, Clone, Eq, PartialEq)]
struct FQP {
    coeffs: Vec<Fq>,
    modulus_coeffs: Vec<U256>,
    degree: usize,
}

impl FQP {
    fn new(coeffs: Vec<Fq>, modulus_coeffs: Vec<U256>) -> Self {
        assert_eq!(coeffs.len(), modulus_coeffs.len());
        Self {
            coeffs: coeffs,
            modulus_coeffs: modulus_coeffs.to_vec(),
            degree: modulus_coeffs.len(),
        }
    }

    fn fq2(coeffs: &[i32]) -> Self {
        assert_eq!(2, coeffs.len());
        Self {
            coeffs: coeffs.iter().map(|c| Fq::new(U256::from(*c))).collect(),
            modulus_coeffs: vec![U256::one(), U256::zero()],
            degree: 2,
        }
    }

    fn fq12(coeffs: &[i32]) -> Self {
        assert_eq!(12, coeffs.len());
        Self {
            coeffs: coeffs.iter().map(|c| Fq::new(U256::from(*c))).collect(),
            modulus_coeffs: fq12_modulus_coeffs(),
            degree: 12,
        }
    }

    fn fq12_one() -> Self {
        let mut v = vec![Fq::new(U256::one())];
        v.resize(12, Fq::new(U256::zero()));
        FQP::fq12_fq(&v)
    }

    fn fq12_fq(coeffs: &[Fq]) -> Self {
        Self {
            coeffs: coeffs.to_vec(),
            modulus_coeffs: fq12_modulus_coeffs(),
            degree: 12,
        }
    }

    fn fq2_g2(p2: [Fq; 2]) -> Self {
        Self {
            coeffs: p2.to_vec(),
            modulus_coeffs: vec![U256::one(), U256::zero()],
            degree: 2,
        }
    }

    fn add(&self, x: &FQP) -> FQP {
        let c = self.coeffs.iter().zip(x.coeffs.iter()).map(|(f1, f2)| f1.add(*f2)).collect();
        FQP {
            coeffs: c,
            modulus_coeffs: self.modulus_coeffs.clone(),
            degree: self.degree,
        }
    }

    fn sub(&self, x: &FQP) -> FQP {
        let c = self.coeffs.iter().zip(x.coeffs.iter()).map(|(f1, f2)| f1.sub(*f2)).collect();
        FQP {
            coeffs: c,
            modulus_coeffs: self.modulus_coeffs.clone(),
            degree: self.degree,
        }
    }

    fn div_fq(&self, x: Fq) -> FQP {
        let ncoeffs: Vec<Fq> = self.coeffs.iter().map(|c| c.div(x)).collect();
        FQP {
            coeffs: ncoeffs,
            modulus_coeffs: self.modulus_coeffs.clone(),
            degree: self.degree,
        }
    }

    fn div_fqp(&self, x: FQP) -> FQP {
        assert_eq!(self.degree, x.degree);
        let v = x.inv();
        let res = self.mul_fqp(&v);
        res
    }

    fn inv(&self) -> FQP {
        fn deg(p: &Vec<Fq>) -> usize {
            let mut d = p.len() - 1;
            while p[d].0.is_zero() && d > 0 {
                d -= 1;
            }
            d
        }

        fn poly_rounded_div(a: &Vec<Fq>, b: &Vec<Fq>) -> Vec<Fq> {
            let dega = deg(&a);
            let degb = deg(&b);
            let mut temp = a.clone();
            let mut o = vec![Fq::new(U256::zero()); a.len()];
            let dd = (dega as i32) - (degb as i32) + 1;
            assert!(dd >= 0);
            for i in (0..dd).rev() {
                let i = i as usize;
                o[i] = o[i].add(temp[degb + i].div(b[degb]));
                for c in 0..(degb + 1) {
                    temp[c + 1] = temp[c + 1].sub(o[c]);
                }
            }
            assert!(deg(&o) + 1 <= o.len());
            o.resize(deg(&o) + 1, Default::default());
            o
        }

        let mut lm = vec![Fq::new(U256::one())];
        lm.resize(self.degree, Fq::new(U256::zero()));
        let mut hm = vec![Fq::new(U256::zero())];
        hm.resize(self.degree, Fq::new(U256::one()));
        let mut low = self.coeffs.clone();
        low.push(Fq::new(U256::zero()));
        let mut high = self.modulus_coeffs.iter().map(|x| Fq::new(*x)).collect();
        while deg(&low) > 0 {
            let mut r = poly_rounded_div(&high, &low);
            assert!(self.degree + 1 - r.len() >= 0);
            r.resize(self.degree + 1, Fq::new(U256::zero()));

            let mut nm = hm.clone();
            let mut new = high.clone();
            assert_eq!(lm.len(), hm.len());
            assert_eq!(lm.len(), low.len());
            assert_eq!(lm.len(), high.len());
            assert_eq!(lm.len(), nm.len());
            assert_eq!(lm.len(), new.len());
            assert_eq!(lm.len(), self.degree + 1);
            for i in 0..self.degree + 1 {
                for j in 0..(self.degree + 1 - i) {
                    nm[i + j] = nm[i + j].sub(lm[i].mul(r[j]));
                    new[i + j] = new[i + j].sub(low[i].mul(r[j]));
                }
            }
            hm = lm;
            lm = nm;
            high = low;
            low = new;
        }
        self.div_fq(low[0])
    }

    fn mul_fq(&self, x: Fq) -> FQP {
        FQP {
            coeffs: self.coeffs.iter().map(|c| c.mul(x)).collect(),
            modulus_coeffs: self.modulus_coeffs.clone(),
            degree: self.degree,
        }
    }

    fn mul_fqp(&self, x: &FQP) -> FQP {
        assert_eq!(self.degree, x.degree);
        let mut b = vec![];
        for i in 0..(self.degree * 2 - 1) {
            b.push(Fq::new(U256::zero()));
        }
        for i in 0..self.degree {
            for j in 0..self.degree {
                b[i + j] = b[i + j].add(self.coeffs[i].mul(x.coeffs[j]));
            }
        }
        while b.len() > self.degree {
            let exp = b.len() - self.degree - 1;
            let top = b.pop().unwrap();
            for i in 0..self.degree {
                b[exp + i] = b[exp + i].sub(top.mul(Fq::new(self.modulus_coeffs[i])));
            }
        }
        let res = FQP { coeffs: b, modulus_coeffs: self.modulus_coeffs.clone(), degree: self.degree };
        res
    }

    fn one(&self) -> Self {
        let mut v = vec![Fq::new(U256::one())];
        v.resize(self.degree, Fq::new(U256::zero()));
        Self {
            coeffs: v,
            modulus_coeffs: self.modulus_coeffs.clone(),
            degree: self.degree,
        }
    }

    fn pow_biguint(&self, p: &BigUint) -> FQP {
        if p.is_zero() {
            let mut v = vec![Fq::new(U256::one())];
            v.resize(self.degree, Fq::new(U256::zero()));
            FQP { coeffs: v, modulus_coeffs: self.modulus_coeffs.clone(), degree: self.degree }
        } else if p.is_one() {
            FQP { coeffs: self.coeffs.clone(), modulus_coeffs: self.modulus_coeffs.clone(), degree: self.degree }
        } else if p.is_even() {
            let p2 : BigUint = p >> 1;
            self.mul_fqp(&self).pow_biguint(&p2)
        } else {
            let p2 : BigUint = p >> 1;
            self.mul_fqp(&self).pow_biguint(&p2).mul_fqp(&self)
        }
    }

    fn pow(&self, p: U256) -> FQP {
        if p.is_zero() {
            let mut v = vec![Fq::new(U256::one())];
            v.resize(self.degree, Fq::new(U256::zero()));
            FQP { coeffs: v, modulus_coeffs: self.modulus_coeffs.clone(), degree: self.degree }
        } else if p == U256::one() {
            FQP { coeffs: self.coeffs.clone(), modulus_coeffs: self.modulus_coeffs.clone(), degree: self.degree }
        } else if p.bit(0) {
            self.mul_fqp(self).pow(p >> 1)
        } else {
            self.mul_fqp(self).pow(p >> 1).mul_fqp(self)
        }
    }

    fn neg(&self) -> FQP {
        let mut v = vec![Fq::new(U256::zero());self.degree];
        for i in 0..self.degree {
            v[i] = self.coeffs[i].neg();
        }
        FQP { coeffs: v, modulus_coeffs: self.modulus_coeffs.clone(), degree: self.degree }
    }
}

fn main() {
    /*
    let a = U256::from_dec_str("13110391464550333261117142677974698181990097963069030438958550202690327139681").unwrap();
    let b = U256::from_dec_str("21888242871839275222246405745257275088696311157297823662689037894645226208583").unwrap();
    let mut lm = U256WithSign { x: U256::one(), is_neg: false };
    let mut hm = U256WithSign { x: U256::zero(), is_neg: false };
    let mut low = U256WithSign { x: a % b, is_neg: false };
    let mut high = U256WithSign { x: b, is_neg: false };
    println!("inv:: a: {}, n: {}, lm: {}, hm: {}, low: {}, high: {}", a, b, lm, hm, low, high);
    loop {
        let c = low.sub(U256WithSign { x: U256::one(), is_neg: false });
        if c.is_neg || c.x.is_zero() {
            break;
        }
        let r = high.div(low);
        // println!("lm: {:?}, r: {:?}, hm: {:?}",lm,r,hm);
        let q1 = hm;
        let q2 = lm.mul(r);
        let nm = hm.sub(lm.mul(r));
        println!("q1: {}, q2: {}, nm: {}",q1,q2,nm);
        let new = high.sub(low.mul(r));

        high = low;
        low = new;
        hm = lm;
        lm = nm;
        println!("inv:: lm: {}, hm: {}, low: {}, high: {}", lm, hm, low, high);
        assert_eq!(format!("{}",lm),"-1");
        break;
    }
     */

    /*
    x=18107605096562723854700666435595314810062186832481469946731400411002744132475
    y=17499317168194804241681774211615986635343204560183427050823794048667776674132
    bad_res_x = 53730168118548197732355161078755349994613603296287790688142546092
    bad_res_y = 9870421536371311977626828777768261154098086310369995172064536975738326893702
    good_res_x = 5214787552945571192215236478449791155895720709003223478753159186543725682906
    good_res_y = 12553987566452835428544018195402269243789459438178333018500397177432382352376

    pt = (FQ(x),FQ(y))
    print(pt)
    ptmul = multiply(pt,32)
    print(ptmul)
    assert(ptmul[0].n==good_res_x)
    assert(ptmul[1].n==good_res_y)
     */
    /*
    let pt = G1Point {
        x: Fq::new(U256::from_dec_str("18107605096562723854700666435595314810062186832481469946731400411002744132475").unwrap()),
        y: Fq::new(U256::from_dec_str("17499317168194804241681774211615986635343204560183427050823794048667776674132").unwrap()),
    };
    let ptmul = pt.mul_fq(Fq::new(U256::from(33)));
     */
    /*
    assert_eq!(ptmul, G1Point {
        x: Fq::new(U256::from_dec_str("5214787552945571192215236478449791155895720709003223478753159186543725682906").unwrap()),
        y: Fq::new(U256::from_dec_str("12553987566452835428544018195402269243789459438178333018500397177432382352376").unwrap()),
    });
     */
    /*
    assert_eq!(ptmul, G1Point {
        x: Fq::new(U256::from_dec_str("10263488358110025510087226086053915515231863465131356542275859284804201740336").unwrap()),
        y: Fq::new(U256::from_dec_str("6902272437764962319374951619802039911156973454275339125186672208635652758702").unwrap()),
    });
    println!("hello");
     */
    /*
    let pt = G1Point {
        x: Fq::new(U256::from_dec_str("18107605096562723854700666435595314810062186832481469946731400411002744132475").unwrap()),
        y: Fq::new(U256::from_dec_str("17499317168194804241681774211615986635343204560183427050823794048667776674132").unwrap()),
    };
    let dpt = pt.double();
    assert_eq!(dpt, G1Point {
        x: Fq::new(U256::from_dec_str("5442964062532895299701493861489680630545642782552326115080478423671196018235").unwrap()),
        y: Fq::new(U256::from_dec_str("14834934525227908230884094748597268296879270652510348750321292389097402766194").unwrap()),
    });
    println!("hello");
     */
    let a1 = U256::from_dec_str("14381381155403816334798493340013261768598472201420748168820652964797987547462").unwrap();
    let a2 = U256::from_dec_str("3018836768243539726548019185072612787535488372041565734163560895787132413963").unwrap();
    let b11 = U256::from_dec_str("1222372640101605630716805870247697620367955331142053227617023362218224623675").unwrap();
    let b12 = U256::from_dec_str("20145010989865280903946427053568635653229019267312571214695898029471466643189").unwrap();
    let b21 = U256::from_dec_str("18627856813852125933293899443711247288107709517554840445817048145846492416046").unwrap();
    let b22 = U256::from_dec_str("7179793276317950049211890536373035911436265064183990876740337674144380042597").unwrap();
    let c1 = U256::from_dec_str("12411515119042684284347439080261290708250768953338792298174140687146844441182").unwrap();
    let c2 = U256::from_dec_str("20138364572167484814256846272295865295005456643842414091875980986217735719454").unwrap();
    let a = [Fq::new(a1), Fq::new(a2)];
    let b = [[Fq::new(b11), Fq::new(b12)], [Fq::new(b21), Fq::new(b22)]];
    let c = [Fq::new(c1), Fq::new(c2)];
    let input = [U256::from_dec_str("32").unwrap()];
    let res = Verifier::verifyProof(a, b, c, &input);
    assert!(!res);
}

/*

self: G1Point(
        Fq(19243360343105341436538200894543113894037152317061756822922778657861755751094),
        Fq(14018277940397341683014465763091036605547806946341036878573379237255260306037))
p2:   G1Point(
        Fq(18107605096562723854700666435595314810062186832481469946731400411002744132475),
        Fq(17499317168194804241681774211615986635343204560183427050823794048667776674132))
l:     Fq(740504438677335617567398191881109606317777637674723293810329998632767455936)
new_x: Fq(20615349946987316985872226152635538518122223301918292800482049448750485629093)
new_y: Fq(7360325878448591122563611461263507904920740457150838477929648565322871406219)
*/
