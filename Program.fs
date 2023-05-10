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

let search (left_point : Point) (right_point : Point) foo m M =
    let x: double = 0.5 * (right_point.X + left_point.X) - (right_point.Z - left_point.Z) / (2.0 * m)
    let z: double = foo x
    let new_point: Point = Point(x, z)
    let prevM: double = abs ((right_point.Z - left_point.Z) / (right_point.X - left_point.X))
    let newM1: double = abs((z - left_point.Z) / (x - left_point.X));
    let newM2: double = abs((right_point.Z - z) / (right_point.X - x));
    if prevM <> M then
        if (M >= max newM1 newM2) then
            let prevR: double = m * (right_point.X - left_point.X) + (right_point.Z - left_point.Z) * (right_point.Z - left_point.Z) / (m * (right_point.X - left_point.X)) - 2.0 * (right_point.Z + left_point.Z);
            let newR1: double = m * (x - left_point.X) + (z - left_point.Z) * (z - left_point.Z) / (m * (x - left_point.X)) - 2.0 * (z + left_point.Z);
            let newR2: double = m * (right_point.X - x) + (right_point.Z - z) * (right_point.Z - z) / (m * (right_point.X - x)) - 2.0 * (right_point.Z + z);
            0, new_point, prevM, newM1, newM2, prevR, newR1, newR2
        else
            2, new_point, prevM, newM1, newM2, 0, 0, 0
    else
        1, new_point, prevM, newM1, newM2, 0, 0, 0



let recalcNewRmap (m: double) (w : List<Point>) =
    let rec loop (m: double) (i: int) (w : List<Point>) (rmap: Map<double,Point>) = 
        if i >= (w.Length - 1) then rmap
        else 
            let i_1_point: Point = w |> List.item i
            let i_point: Point = w |> List.item (i + 1)
            let R : double = m * (i_point.X - i_1_point.X) + (i_point.Z - i_1_point.Z) * (i_point.Z - i_1_point.Z) / (m * (i_point.X - i_1_point.X)) - 2.0 * (i_point.Z + i_1_point.Z)
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
        let newM: double = List.max mset
        let newm: double = calc_m newM r
        recalcNewRmap newm w, newM, newm

let insert_point_parallel (new_point: Point) m M mset (rmap: Map<double, Point>) (w : List<Point>) newM1 newM2 Mprev Rprev newR1 newR2 =
    let w_new = (List.append w [new_point]) |> List.sortBy (fun elem -> elem.X)
    let new_elem_idx = w_new |> List.findIndex (fun x -> x.X.Equals new_point.X)
    let new_elem = w_new |> List.item new_elem_idx
    let after_new_elem = w_new |> List.item (new_elem_idx + 1)
    let Mset_new : List<double> = 
        (mset |> List.filter (fun (x: double) -> x <> Mprev)) 
        |> List.append [newM1; newM2]
    if (M = List.max Mset_new) then 
        let new_rmap: Map<double, Point>  = rmap |> Map.remove Rprev |> Map.add newR1 new_elem |> Map.add newR2 after_new_elem
        w_new, Mset_new, new_rmap
    else 
        w_new, Mset_new, rmap

let new_interval (rmap: Map<double, Point>) (w : List<Point>) =
    let max_R: double  = (rmap |> Map.toList) |> List.map fst |> List.max
    let new_t = rmap |> Map.find max_R
    let idx = w |> List.findIndex (fun (x) -> (x : Point).X = new_t.X)
    let new_t_1 = w |> List.item (idx - 1)
    new PairOfPoints(new_t, new_t_1)

let mainFunction (PairOfPoints: PairOfPoints) (m: double) M mset (rmap: Map<double, Point>) (w : List<Point>) res r foo =
    let x: double = 0.5 * (PairOfPoints.T.X + PairOfPoints.T_1.X) -  (PairOfPoints.T.Z - PairOfPoints.T_1.Z) / (2.0 * m)
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
    let max_R: double  = (new_rmap |> Map.toList) |> List.map fst |> List.max
    let new_t = new_rmap |> Map.find max_R
    let idx = w_new |> List.findIndex (fun (x) -> (x : Point).X = new_t.X)
    let new_t_1 = w_new |> List.item (idx - 1)
    new PairOfPoints(new_t, new_t_1), w_new, newM, newm, new_rmap, Mset_new, new_res

let pickNMax (rmap: Map<double, Point>) (n: int)=
    let rec loop  (rmap: Map<double, Point>) (i: int) (rmax : List<double * Point>) = 
        if (i = n) then
            rmax
        else
            let max_R: double  = (rmap |> Map.toList) |> List.map fst |> List.max
            let point = rmap |> Map.find max_R
            let new_rmax: List<double * Point> =  rmax @ [(max_R, point)]
            loop (rmap |> Map.remove max_R) (i+1) new_rmax
    loop rmap 0 List.Empty

let rec iter (PairOfPoints: PairOfPoints) epsilon n nmax m M mset rmap (w: List<Point>) res r foo = 

    let nthread = 4
    if abs(PairOfPoints.T.X - PairOfPoints.T_1.X) <= epsilon or n > nmax then
        n, res
    else if n < 4 then
        let (PairOfPoints2: PairOfPoints, new_w, newM, new_m, new_rmap, new_mset, new_res) = mainFunction PairOfPoints m M mset rmap w res r foo
        iter PairOfPoints2 epsilon (n+1) nmax new_m newM new_mset new_rmap new_w new_res r foo
        else 
            let Rmax : List<double * Point> = pickNMax rmap nthread
            let return_codes: int array = Array.zeroCreate nthread
            let new_point_arr: Point array = Array.zeroCreate nthread
            let prevM_arr: double array = Array.zeroCreate nthread
            let newM1_arr: double array = Array.zeroCreate nthread
            let newM2_arr: double array = Array.zeroCreate nthread
            let prevR_arr: double array = Array.zeroCreate nthread
            let newR1_arr: double array = Array.zeroCreate nthread
            let newR2_arr: double array = Array.zeroCreate nthread
            Parallel.For(0, nthread, fun (i: int) ->
                let (key: double, value) =  Rmax.[i]
                let idx = List.findIndex (fun (x: Point)-> x.X = value.X) w
                let (return_code, new_point, prevM, newM1, newM2, prevR, newR1, newR2) = search (List.item (idx-1) w)  (List.item idx w) foo m M
                return_codes.[i] <- return_code
                prevM_arr.[i] <- prevM
                newM1_arr.[i] <- newM1
                newM2_arr.[i] <- newM2
                prevR_arr.[i] <- prevR
                newR1_arr.[i] <- newR1
                newR2_arr.[i] <- newR2
                new_point_arr.[i] <- new_point
            )
            if (Array.sum return_codes) > 0 then 
                let i = Array.findIndex (fun x -> x > 0) return_codes 
                let (new_w, new_mset, dmp) = 
                    insert_point_parallel new_point_arr.[i] m M mset rmap w newM1_arr.[i] newM2_arr.[i] prevM_arr.[i] prevR_arr.[i] newR1_arr.[i] newR2_arr.[i]   
   
                let newM = List.max new_mset
                let new_m = calc_m newM r
                let new_rmap = recalcNewRmap new_m new_w
                let new_res = minOfPrevAndNewRes res new_point_arr.[i]
                iter (new_interval new_rmap new_w) epsilon (n+1) nmax new_m newM new_mset new_rmap new_w new_res r foo
            else 
                let rec loop nthread m M (i: int) (w : List<Point>) (rmap: Map<double,Point>) mset res (newM1_arr : double array) (newM2_arr : double array) (prevM_arr : double array) (prevR_arr : double array) (newR1_arr : double array) (newR2_arr : double array) = 
                    if i >= nthread then w, mset, rmap, res
                    else 
                        let (new_w, new_mset, new_rmap) = 
                            insert_point_parallel new_point_arr.[i] m M mset rmap w newM1_arr.[i] newM2_arr.[i] prevM_arr.[i] prevR_arr.[i] newR1_arr.[i] newR2_arr.[i]
                        let new_res = minOfPrevAndNewRes res new_point_arr.[i]
                        loop nthread m M (i+1) new_w new_rmap new_mset new_res newM1_arr newM2_arr prevM_arr prevR_arr newR1_arr newR2_arr
                let (new_w, new_mset, new_rmap, new_res) = 
                    loop nthread m M 0 w rmap mset res newM1_arr newM2_arr prevM_arr prevR_arr newR1_arr newR2_arr
                
                iter (new_interval new_rmap new_w) epsilon (n+1) nmax m M new_mset new_rmap new_w new_res r foo


let GSA_parallel_v2 foo (a: double) (b: double) (r: double) (epsilon: double) (nmax: int) =
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
        let (n: int, res : Point) = GSA_parallel_v2 (ArrayFunc.Item(i)) (bounds.Item(j)) (bounds.Item(j+1)) 2.5 0.01 10000
        timer.Stop()
        printfn "n = %O" n
        printfn "elapsed = %O" timer.Elapsed
        printfn "res.X  = %f, res.Z = %f" res.X res.Z
        exec_loop (i+1) (j+2)

let ret = exec_loop 0 0
