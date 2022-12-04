open System

type Point() =
   let mutable x = 0.0
   let mutable z = 0.0
   member this.X
      with get() = x
      and set(value) = x <- value
   member this.Z
      with get() = z
      and set(value) = z <- value
   interface IComparable with
     member this.CompareTo(x2: obj) =
       if this.X <= (x2 :?> Point).X 
         then 1
         else 0 


type Pair() =
   let mutable t = new Point()
   let mutable t_1 = new Point()
   member this.T
      with get() = t
      and set(value) = t <- value
   member this.T_1
      with get() = t_1
      and set(value) = t_1 <- value

let NUM_ITER = 1000000.0
let M_PI = 3.14159265358979323846
let sincos x = Math.Sin(x)*Math.Sin(x) + Math.Cos(x)*Math.Cos(x)

let calculate items =
    items |> List.map sincos

let calc_m (M:float) (r:float) =
    if M > 0.0 then r*M
    else 1.0

let find_res (prev : Point) (new_p : Point) =
    if new_p.Z < prev.Z then new_p
    else prev

let calc_new_Rmap m (w : List<Point>) =
    let rec loop m i (w : List<Point>) rmap = 
        if i >= (w.Length - 1) then rmap
        else 
            let i_1_point = w |> List.item i
            let i_point = w |> List.item (i + 1)
            let R = m * (i_point.X - i_1_point.X) + (i_point.Z - i_1_point.Z) * (i_point.Z - i_1_point.Z) / (m * (i_point.X - i_1_point.X)) - 2.0 * (i_point.Z + i_1_point.Z)
            let new_rmap : Map<float, Point> = rmap |> Map.add R i_point
            loop m (i + 1) w new_rmap
    loop m 0 w Map.empty

let recalc_rmap mset (rmap : Map<float, Point>) m M (new_elem: Point) (before_new_elem: Point) (after_new_elem: Point) (w : List<Point>) r =
    if (M = List.max mset) then 
        let dmp = rmap |> Map.remove (m * (after_new_elem.X - before_new_elem.X) + 
        (after_new_elem.Z - before_new_elem.Z) * (after_new_elem.Z - before_new_elem.Z) / 
        (m * (after_new_elem.X - before_new_elem.X)) - 2.0 * (after_new_elem.Z + before_new_elem.Z))
        let dmp_R = m * (after_new_elem.X - new_elem.X) + (after_new_elem.Z - new_elem.Z) * (after_new_elem.Z - new_elem.Z) / (m * (after_new_elem.X - new_elem.X)) - 2.0 * (after_new_elem.Z + new_elem.Z)
        let dmp_R2 = m * (new_elem.X - before_new_elem.X) + (new_elem.Z - before_new_elem.Z) * (new_elem.Z - before_new_elem.Z) / (m * (new_elem.X - before_new_elem.X)) - 2.0 * (new_elem.Z + before_new_elem.Z)
        let dmp3 = dmp |> Map.add dmp_R after_new_elem
        dmp3 |> Map.add dmp_R2 new_elem, M, m
    else 
        let newM = List.max mset
        let newm = calc_m newM r
        calc_new_Rmap newm w, newM, newm

let new_points (pair: Pair) m M mset (rmap: Map<float, Point>) (w : List<Point>) res r foo =
    let x = 0.5 * (pair.T.X + pair.T_1.X) -  (pair.T.Z - pair.T_1.Z) / (2.0 * m)
    let new_point = new Point(X = x, Z = foo x)
    let new_res = find_res res new_point
    let w_new = (List.append w [new_point]) |> List.sortBy (fun elem -> elem.X)
    let new_elem_idx = w_new |> List.findIndex (fun x -> x.X.Equals new_point.X) 
    let new_elem = w_new |> List.item new_elem_idx
    let before_new_elem = w_new |> List.item (new_elem_idx - 1)
    let after_new_elem = w_new |> List.item (new_elem_idx + 1)
    let Mset_new : List<float> = 
        (mset |> List.filter (fun x -> x <> abs ((after_new_elem.Z - before_new_elem.Z) / (after_new_elem.X - before_new_elem.X)))) 
        |> List.append [abs ((after_new_elem.Z - new_elem.Z) / (after_new_elem.X - new_elem.X));
        abs ((new_elem.Z - before_new_elem.Z) / (new_elem.X - before_new_elem.X))]
    let (new_rmap : Map<float, Point>, newM, newm) = recalc_rmap Mset_new rmap m M new_elem before_new_elem after_new_elem w_new r
    let max_R  = (new_rmap |> Map.toList) |> List.map fst |> List.max
    let new_t = new_rmap |> Map.find max_R
    let idx = w_new |> List.findIndex (fun (x) -> (x : Point).X = new_t.X)
    let new_t_1 = w_new |> List.item (idx - 1)
    new Pair(T = new_t, T_1 = new_t_1), w_new, newM, newm, new_rmap, Mset_new, new_res


let rec iter (pair: Pair) epsilon n nmax m M mset rmap w res r foo = 
    if pair.T.X - pair.T_1.X <= epsilon then res
    else if n > nmax then 
        printf "n = %d" n
        res
    else 
        let (pair2: Pair, new_w, newM, new_m, new_rmap, new_mset, new_res) = new_points pair m M mset rmap w res r foo
        iter pair2 epsilon (n+1) nmax new_m newM new_mset new_rmap new_w new_res r foo



let GSA foo a b r epsilon nmax =
    let point_a = new Point(X = a, Z = foo a)
    let point_b = new Point(X = b, Z = foo b)
    let w = [point_a; point_b]   
    let M = abs ((foo b - foo a) / (b - a))
    let Mset = [M]
    let m = calc_m M r
    let R = m * (b - a) + (foo b - foo a) * (foo b - foo a) / (m * (b - a)) - 2.0 * (foo b + foo a)
    let Rset = Map.empty.Add(R, point_b)
    let res = find_res point_a point_b
    iter (new Pair(T = point_b, T_1 = point_a)) epsilon 0 nmax m M Mset Rset w res r foo


let foo0 x = sin x + sin (10.0 * x / 3.0) // 2.7 7.5 // res.X  = 5.148005, res.Z = -1.899569
let foo1 x = 2.0 * (x - 3.0) * (x - 3.0) + exp (x * x / 2.0) // -3.0 3.0 // res.X  = 1.589118, res.Z = 7.515945
let foo2 x = -1.0 * (1.0 * sin ((1.0 + 1.0) * x + 1.0) + 
    2.0 * sin ((2.0 + 1.0) * x + 2.0) + 
    3.0 * sin ((3.0 + 1.0) * x + 3.0) +
    4.0 * sin ((4.0 + 1.0) * x + 4.0) + 
    5.0 * sin ((5.0 + 1.0) * x + 5.0))  // 0.0 10.0 // res.X  = 5.793274, res.Z = -12.030911
let foo3 x = (3.0 * x - 1.4) * sin (18.0 * x) // 0.0 1.2 // res.X  = 0.967300, res.Z = -1.488708
let foo4 x = -1.0 * (x + sin x) * exp (-1.0 * x * x) // -10.0 10.0 // res.X  = 0.682491, res.Z = -0.824224
let foo5 x = sin x + sin (10.0 * x / 3.0) + log x - 0.84 * x + 3.0 // 2.7 7.5 // res.X  = 5.202703, res.Z = -1.601256
let foo6 x = -1.0 * sin (2.0 * M_PI * x) * exp (-1.0 * x) // 0.0 4.0 // res.X  = 0.226858, res.Z = -0.788623
let foo7 x = (x * x - 5.0 * x + 6.0) / (x * x + 1.0) // -5.0 5.0 // res.X  = 2.413894, res.Z = -0.035534
let foo8 x = -1.0 * x + sin (3.0 * x) - 1.0  // 0.0 6.5 // res.X  = 5.876861, res.Z = -7.815607

let res = GSA foo8 0.0 6.5 2.0 0.01 1000
printfn "res.X  = %f, res.Z = %f" res.X res.Z