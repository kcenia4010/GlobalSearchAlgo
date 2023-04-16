open System
open System.Threading.Tasks

type Point(x:double, z:double) =
  member this.X = x
  member this.Z = z
  interface IComparable with
     member this.CompareTo(x2: obj) =
       if this.X <= (x2 :?> Point).X
         then 1
         else 0

type PairOfPoints(t:Point, t_1:Point) =
   member this.T = t
   member this.T_1 = t_1

let NUM_ITER = 1000000.0
let M_PI =3.14159265358979323846264338327950288 
let sincos (x: double) = Math.Sin(x)*Math.Sin(x) + Math.Cos(x)*Math.Cos(x)

let rec slow_down (item: double) (i: double) =
    if i >= NUM_ITER then item
    else slow_down (item * (sincos i)) (i + 1.0)

let calc_m (M:double) (r:double) =
    if M > 0.0 then r*M
    else 1.0

let minOfPrevAndNewRes (prev : Point) (new_p : Point) =
    if new_p.Z < prev.Z then new_p
    else prev

let recalcNewRmap (m: double) (w : List<Point>) =
    let rec loop (m: double) (i: int) (w : List<Point>) (rmap: Map<double,Point>) = 
        if i >= (w.Length - 1) then rmap
        else 
            let i_1_point: Point = w |> List.item i
            let i_point: Point = w |> List.item (i + 1)
            let R = m * (i_point.X - i_1_point.X) + (i_point.Z - i_1_point.Z) * (i_point.Z - i_1_point.Z) / (m * (i_point.X - i_1_point.X)) - 2.0 * (i_point.Z + i_1_point.Z)
            let new_rmap : Map<double, Point> = rmap |> Map.add R i_point
            loop m (i + 1) w new_rmap
    loop m 0 w Map.empty

let calcRmap mset (rmap : Map<double, Point>) m M (new_elem: Point) (before_new_elem: Point) (after_new_elem: Point) (w : List<Point>) r =
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
        recalcNewRmap newm w, newM, newm

let mainFunction (PairOfPoints: PairOfPoints) m M mset (rmap: Map<double, Point>) (w : List<Point>) res r foo =
    let x = 0.5 * (PairOfPoints.T.X + PairOfPoints.T_1.X) -  (PairOfPoints.T.Z - PairOfPoints.T_1.Z) / (2.0 * m)
    let new_point = new Point(x, foo x)
    let new_res = minOfPrevAndNewRes res new_point
    let w_new = (List.append w [new_point]) |> List.sortBy (fun elem -> elem.X)
    let new_elem_idx = w_new |> List.findIndex (fun x -> x.X.Equals new_point.X) 
    let new_elem = w_new |> List.item new_elem_idx
    let before_new_elem = w_new |> List.item (new_elem_idx - 1)
    let after_new_elem = w_new |> List.item (new_elem_idx + 1)
    let Mset_new : List<double> = 
        (mset |> List.filter (fun x -> x <> abs ((after_new_elem.Z - before_new_elem.Z) / (after_new_elem.X - before_new_elem.X)))) 
        |> List.append [abs ((after_new_elem.Z - new_elem.Z) / (after_new_elem.X - new_elem.X));
        abs ((new_elem.Z - before_new_elem.Z) / (new_elem.X - before_new_elem.X))]
    let (new_rmap : Map<double, Point>, newM, newm) = calcRmap Mset_new rmap m M new_elem before_new_elem after_new_elem w_new r
    let max_R  = (new_rmap |> Map.toList) |> List.map fst |> List.max
    let new_t = new_rmap |> Map.find max_R
    let idx = w_new |> List.findIndex (fun (x) -> (x : Point).X = new_t.X)
    let new_t_1 = w_new |> List.item (idx - 1)
    new PairOfPoints(new_t, new_t_1), w_new, newM, newm, new_rmap, Mset_new, new_res


let rec iter (PairOfPoints: PairOfPoints) epsilon n nmax m M mset rmap w res r foo = 
    if PairOfPoints.T.X - PairOfPoints.T_1.X <= epsilon or n > nmax then 
        n, res
    else 
        let (PairOfPoints2: PairOfPoints, new_w, newM, new_m, new_rmap, new_mset, new_res) = mainFunction PairOfPoints m M mset rmap w res r foo
        iter PairOfPoints2 epsilon (n+1) nmax new_m newM new_mset new_rmap new_w new_res r foo

let GSA foo a b r epsilon nmax =
    let point_a = new Point(a, foo a)
    let point_b = new Point(b, foo b)
    let w = [point_a; point_b]   
    let M = abs ((foo b - foo a) / (b - a))
    let Mset = [M]
    let m = calc_m M r
    let R = m * (b - a) + (foo b - foo a) * (foo b - foo a) / (m * (b - a)) - 2.0 * (foo b + foo a)
    let Rset = Map.empty.Add(R, point_b)
    let res = minOfPrevAndNewRes point_a point_b
    iter (new PairOfPoints(point_b, point_a)) epsilon 0 nmax m M Mset Rset w res r foo

let rec findMin (arr: Point array) (min: Point) (index: int) =
    if index >= arr.Length then
        min
    else
        let current = arr.[index]
        let newMin = if current.Z < min.Z then current else min
        findMin arr newMin (index + 1)

let rec findMax (arr: int array) (max: int) (index: int) =
    if index >= arr.Length then
        max
    else
        let current = arr.[index]
        let newMax = if current > max then current else max
        findMax arr newMax (index + 1)

let GSA_parallel_v1 foo (a:double) (b:double) r epsilon (nmax: int)  =
    let ntrhead = 4
    let dist: double = abs (b - a)
    let result_points: Point array = Array.zeroCreate ntrhead
    let result_n: int array = Array.zeroCreate ntrhead
    Parallel.For(0, ntrhead, fun (i: int) ->
        let start: double = a + double i * dist/double ntrhead
        let finish: double = a + (double i + 1.0) * dist/double ntrhead
        let (n1 : int, res1 : Point) = GSA foo start finish r epsilon (nmax/ntrhead)
        result_points.[i] <- res1
        result_n.[i] <- n1
    )
    findMax (result_n) result_n.[0] 0, findMin (result_points) result_points.[0] 0

let foo0 x = sin x + sin (10.0 * x / 3.0) * slow_down 1.0 0.0// 2.7 7.5 // res.X  = 5.148005, res.Z = -1.899569
let foo1 x = 2.0 * (x - 3.0) * (x - 3.0) + exp (x * x / 2.0) * slow_down 1.0 0.0// -3.0 3.0 // res.X  = 1.589118, res.Z = 7.515945
let foo2 x = -1.0 * (1.0 * sin ((1.0 + 1.0) * x + 1.0) + 
    2.0 * sin ((2.0 + 1.0) * x + 2.0) + 
    3.0 * sin ((3.0 + 1.0) * x + 3.0) +
    4.0 * sin ((4.0 + 1.0) * x + 4.0) + 
    5.0 * sin ((5.0 + 1.0) * x + 5.0)) * slow_down 1.0 0.0  // 0.0 10.0 // res.X  = 5.793274, res.Z = -12.030911
let foo3 x = (3.0 * x - 1.4) * sin (18.0 * x) * slow_down 1.0 0.0 // 0.0 1.2 // res.X  = 0.967300, res.Z = -1.488708
let foo4 x = -1.0 * (x + sin x) * exp (-1.0 * x * x) * slow_down 1.0 0.0 // -10.0 10.0 // res.X  = 0.682491, res.Z = -0.824224
let foo5 x = sin x + sin (10.0 * x / 3.0) + log x - 0.84 * x + 3.0 * slow_down 1.0 0.0 // 2.7 7.5 // res.X  = 5.202703, res.Z = -1.601256
let foo6 x = -1.0 * sin (2.0 * M_PI * x) * exp (-1.0 * x) * slow_down 1.0 0.0 // 0.0 4.0 // res.X  = 0.226858, res.Z = -0.788623
let foo7 x = (x * x - 5.0 * x + 6.0) / (x * x + 1.0) * slow_down 1.0 0.0 // -5.0 5.0 // res.X  = 2.413894, res.Z = -0.035534
let foo8 x = -1.0 * x + sin (3.0 * x) - 1.0 * slow_down 1.0 0.0// 0.0 6.5 // res.X  = 5.876861, res.Z = -7.815607

let ArrayFunc = [foo0; foo1; foo2; foo3; foo4; foo5; foo6; foo7; foo8]
let bounds = [2.7; 7.5; -3.0; 3.0; 0.0; 10.0; 0.0; 1.2; -10.0; 10.0; 2.7; 7.5; 0.0; 4.0; -5.0; 5.0; 0.0; 6.5]

let rec exec_loop i j =
    if i > 8 then 0 
    else
        printfn "function%O" i
        let timer = System.Diagnostics.Stopwatch.StartNew()
        let (n: int, res : Point) = GSA_parallel_v1 (ArrayFunc.Item(i)) (bounds.Item(j)) (bounds.Item(j+1)) 2.5 0.01 10000
        timer.Stop()
        printfn "n = %O" n
        printfn "elapsed = %O" timer.Elapsed
        printfn "res.X  = %f, res.Z = %f" res.X res.Z
        exec_loop (i+1) (j+2)

let ret = exec_loop 0 0
