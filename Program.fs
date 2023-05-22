open System
open System.Threading.Tasks
open System.Runtime.InteropServices

[<DllImport("grishagin.dll", CallingConvention = CallingConvention.Cdecl)>]
extern void CreateGrishaginFunction()

[<DllImport("grishagin.dll", CallingConvention = CallingConvention.Cdecl)>]
extern void MySetFunctionNumber(int i)

[<DllImport("grishagin.dll", CallingConvention = CallingConvention.Cdecl)>]
extern double GetBoundsA0()

[<DllImport("grishagin.dll", CallingConvention = CallingConvention.Cdecl)>]
extern double GetBoundsA1() 

[<DllImport("grishagin.dll", CallingConvention = CallingConvention.Cdecl)>]
extern double GetBoundsB0()

[<DllImport("grishagin.dll", CallingConvention = CallingConvention.Cdecl)>]
extern double GetBoundsB1() 

[<DllImport("grishagin.dll", CallingConvention = CallingConvention.Cdecl)>]
extern double MyCalculate(double y0, double y1)

[<DllImport("grishagin.dll", CallingConvention = CallingConvention.Cdecl)>]
extern void FreeGrishaginFunction()



[<DllImport("evolvent.dll", CallingConvention = CallingConvention.Cdecl)>]
extern void CreateEvovlent(int idx)

[<DllImport("evolvent.dll", CallingConvention = CallingConvention.Cdecl)>]
extern void CopyEvovlent(int idx)

[<DllImport("evolvent.dll", CallingConvention = CallingConvention.Cdecl)>]
extern void SetBounds(int idx, double A0, double A1, double B0, double B1)

[<DllImport("evolvent.dll", CallingConvention = CallingConvention.Cdecl)>]
extern double GetImage0(int idx, double x)

[<DllImport("evolvent.dll", CallingConvention = CallingConvention.Cdecl)>]
extern double GetImage1(int idx, double x)

[<DllImport("evolvent.dll", CallingConvention = CallingConvention.Cdecl)>]
extern void FreeEvolvent(int idx)

type Point(x:double, z:double) =
    member this.X = x
    member this.Z = z
    override this.Equals(obj: obj) =
        match obj with
        | :? Point as other -> this.X = other.X && this.Z = other.Z
        | _ -> false
    override this.GetHashCode() =
        hash (this.X, this.Z)
    interface System.IComparable with
        member this.CompareTo(obj : obj) =
            match obj with
            | :? Point as other ->
                if this.X < other.X then -1
                elif this.X = other.X && this.Z = other.Z then 0
                else 1
            | _ -> invalidArg "obj" "The argument must be of type Point"


type PairOfPoints(t:Point, t_1:Point) =
   member this.T = t
   member this.T_1 = t_1

type Result(f:double, y:List<double>) =
    member this.Y = y
    member this.F = f

let NUM_ITER = 1000000.0
let M_PI =3.14159265358979323846264338327950288 
let sincos (x: double) = Math.Sin(x)*Math.Sin(x) + Math.Cos(x)*Math.Cos(x)

let rec slow_down (item: double) (i: double) =
    if i >= NUM_ITER then item
    else slow_down (item * (sincos i)) (i + 1.0)

let calc_m (M:double) (r:double) =
    if M > 0.0 then r*M
    else 1.0

let minOfPrevAndNewRes (x : double) (result : Result)  =
    let new_y = [GetImage0(0, x);GetImage1(0, x)]
    let z = MyCalculate((new_y.Item(0)), (new_y.Item(1)))
    if z < result.F then (new Result(z, new_y))
    else result

let search  idx (left_point : Point) (right_point : Point) m M =
    let x: double = 0.5 * (right_point.X + left_point.X) - (right_point.Z - left_point.Z) / (2.0 * m)
    let new_y = [GetImage0(idx, x);GetImage1(idx, x)]
    let new_point = new Point(x, MyCalculate((new_y.Item(0)), (new_y.Item(1))))
    let prevM: double = abs ((right_point.Z - left_point.Z) / (right_point.X - left_point.X))
    let newM1: double = abs((new_point.Z - left_point.Z) / (x - left_point.X));
    let newM2: double = abs((right_point.Z - new_point.Z ) / (right_point.X - x));
    if prevM <> M then
        if (M >= max newM1 newM2) then
            let prevR: double = m * (right_point.X - left_point.X) + (right_point.Z - left_point.Z) * (right_point.Z - left_point.Z) / (m * (right_point.X - left_point.X)) - 2.0 * (right_point.Z + left_point.Z);
            let newR1: double = m * (x - left_point.X) + (new_point.Z  - left_point.Z) * (new_point.Z  - left_point.Z) / (m * (x - left_point.X)) - 2.0 * (new_point.Z  + left_point.Z);
            let newR2: double = m * (right_point.X - x) + (right_point.Z - new_point.Z ) * (right_point.Z - new_point.Z ) / (m * (right_point.X - x)) - 2.0 * (right_point.Z + new_point.Z );
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
    if (M = Set.maxElement mset) then 
        let dmp = rmap |> Map.remove (m * (after_new_elem.X - before_new_elem.X) + 
         (after_new_elem.Z - before_new_elem.Z) * (after_new_elem.Z - before_new_elem.Z) / 
         (m * (after_new_elem.X - before_new_elem.X)) - 2.0 * (after_new_elem.Z + before_new_elem.Z))
        let dmp_R = m * (after_new_elem.X - new_elem.X) + (after_new_elem.Z - new_elem.Z) * (after_new_elem.Z - new_elem.Z) / (m * (after_new_elem.X - new_elem.X)) - 2.0 * (after_new_elem.Z + new_elem.Z)
        let dmp_R2 = m * (new_elem.X - before_new_elem.X) + (new_elem.Z - before_new_elem.Z) * (new_elem.Z - before_new_elem.Z) / (m * (new_elem.X - before_new_elem.X)) - 2.0 * (new_elem.Z + before_new_elem.Z)
        let dmp3 = dmp |> Map.add dmp_R after_new_elem
        dmp3 |> Map.add dmp_R2 new_elem, M, m
    else 
        let newM: double = Set.maxElement mset
        let newm: double = calc_m newM r
        recalcNewRmap newm w, newM, newm

let insert_point_parallel (new_point: Point) m M mset (rmap: Map<double, Point>) (w : List<Point>) newM1 newM2 Mprev Rprev newR1 newR2 =
    let w_new = (List.append w [new_point]) |> List.sortBy (fun elem -> elem.X)
    let new_elem_idx = w_new |> List.findIndex (fun x -> x.X.Equals new_point.X)
    let new_elem = w_new |> List.item new_elem_idx
    let after_new_elem = w_new |> List.item (new_elem_idx + 1)
    let Mset_new : Set<double> = 
        mset |> Set.remove (Mprev)
        |> Set.add (newM1)
        |> Set.add (newM2)
    if (M = Set.maxElement Mset_new) then 
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

let mainFunction (PairOfPoints: PairOfPoints) m M mset (rmap: Map<double, Point>) (w : List<Point>) (res : Result) r =
    let x = 0.5 * (PairOfPoints.T.X + PairOfPoints.T_1.X) -  (PairOfPoints.T.Z - PairOfPoints.T_1.Z) / (2.0 * m)
    let new_y = [GetImage0(0, x);GetImage1(0, x)]
    let new_point = new Point(x, MyCalculate((new_y.Item(0)), (new_y.Item(1))))
    let new_res : Result = minOfPrevAndNewRes new_point.X res
    let w_new = (List.append w [new_point]) |> List.sortBy (fun elem -> elem.X)
    let new_elem_idx = w_new |> List.findIndex (fun x -> x.X.Equals new_point.X) 
    let new_elem = w_new |> List.item new_elem_idx
    let before_new_elem = w_new |> List.item (new_elem_idx - 1)
    let after_new_elem = w_new |> List.item (new_elem_idx + 1)
    let Mset_new : Set<double> = 
        mset |> Set.remove (abs ((after_new_elem.Z - before_new_elem.Z) / (after_new_elem.X - before_new_elem.X)))
        |> Set.add (abs ((after_new_elem.Z - new_elem.Z) / (after_new_elem.X - new_elem.X)))
        |> Set.add (abs ((new_elem.Z - before_new_elem.Z) / (new_elem.X - before_new_elem.X)))
    let (new_rmap : Map<double, Point>, newM, newm) = calcRmap Mset_new rmap m M new_elem before_new_elem after_new_elem w_new r
    let max_R, new_t = Map.maxKeyValue new_rmap
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

let rec iter (PairOfPoints: PairOfPoints) epsilon n nmax m M mset rmap (w: List<Point>) (res: Result) r = 
    let nthread = 4
    for i in 1..nthread do
        CopyEvovlent(i)
    if abs(PairOfPoints.T.X - PairOfPoints.T_1.X) <= epsilon or n > nmax then
        n, res
    else if n < 4 then
        let (PairOfPoints2: PairOfPoints, new_w, newM, new_m, new_rmap, new_mset, new_res) = mainFunction PairOfPoints m M mset rmap w res r
        iter PairOfPoints2 epsilon (n+1) nmax new_m newM new_mset new_rmap new_w new_res r
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
            let task = Task.Run(fun () -> Parallel.For(0, nthread, fun (i: int) ->
                let (key: double, value) =  Rmax.[i]
                let idx = List.findIndex (fun (x: Point)-> x.X = value.X) w
                let (return_code, new_point, prevM, newM1, newM2, prevR, newR1, newR2) = search i (List.item (idx-1) w)  (List.item idx w) m M
                return_codes.[i] <- return_code
                prevM_arr.[i] <- prevM
                newM1_arr.[i] <- newM1
                newM2_arr.[i] <- newM2
                prevR_arr.[i] <- prevR
                newR1_arr.[i] <- newR1
                newR2_arr.[i] <- newR2
                new_point_arr.[i] <- new_point
            ))
            task.Wait()
            for i in 1..nthread do
                FreeEvolvent(i)
            if (Array.sum return_codes) > 0 then 
                let i = Array.findIndex (fun x -> x > 0) return_codes 
                let (new_w, new_mset, dmp) = 
                    insert_point_parallel new_point_arr.[i] m M mset rmap w newM1_arr.[i] newM2_arr.[i] prevM_arr.[i] prevR_arr.[i] newR1_arr.[i] newR2_arr.[i]   
   
                let newM = Set.maxElement new_mset
                let new_m = calc_m newM r
                let new_rmap = recalcNewRmap new_m new_w
                let new_res : Result = minOfPrevAndNewRes new_point_arr.[i].X res
                iter (new_interval new_rmap new_w) epsilon (n+1) nmax new_m newM new_mset new_rmap new_w new_res r
            else 
                let rec loop nthread m M (i: int) (w : List<Point>) (rmap: Map<double,Point>) mset (res: Result) (newM1_arr : double array) (newM2_arr : double array) (prevM_arr : double array) (prevR_arr : double array) (newR1_arr : double array) (newR2_arr : double array) = 
                    if i >= nthread then w, mset, rmap, res
                    else 
                        let (new_w, new_mset, new_rmap) = 
                            insert_point_parallel new_point_arr.[i] m M mset rmap w newM1_arr.[i] newM2_arr.[i] prevM_arr.[i] prevR_arr.[i] newR1_arr.[i] newR2_arr.[i]
                        let new_res : Result = minOfPrevAndNewRes new_point_arr.[i].X res
                        loop nthread m M (i+1) new_w new_rmap new_mset new_res newM1_arr newM2_arr prevM_arr prevR_arr newR1_arr newR2_arr
                let (new_w, new_mset, new_rmap, new_res) = 
                    loop nthread m M 0 w rmap mset res newM1_arr newM2_arr prevM_arr prevR_arr newR1_arr newR2_arr
                
                iter (new_interval new_rmap new_w) epsilon (n+1) nmax m M new_mset new_rmap new_w new_res r


let GSA_parallel_v2 (Abounds: List<double>) (Bbounds: List<double>) (r: double) (epsilon: double) (nmax: int) =
    CreateEvovlent(0)
    SetBounds(0, (Abounds.Item(0)), (Abounds.Item(1)), (Bbounds.Item(0)), (Bbounds.Item(1)))
    let A = 0.0
    let B = 1.0
    let y1 = [GetImage0(0, A);GetImage1(0, A)]
    let y2 = [GetImage0(0, B);GetImage1(0, B)]
    let point_a = new Point(A, MyCalculate((y1.Item(0)), (y1.Item(1))))
    let point_b = new Point(B, MyCalculate((y2.Item(0)), (y2.Item(1))))
    let w = [point_a; point_b]   
    let M = abs ((point_b.Z - point_a.Z) / (B - A))
    let Mset = Set.empty.Add(M)
    let m = calc_m M r
    let R = m * (B - A) + (point_b.Z - point_a.Z) * (point_b.Z - point_a.Z) / (m * (B - A)) - 2.0 * (point_b.Z + point_a.Z)
    let Rset = Map.empty.Add(R, point_b)
    let res : Result = minOfPrevAndNewRes point_b.X (new Result(point_a.Z, y1))
    let (n, result) = iter (new PairOfPoints(point_b, point_a)) epsilon 0 nmax m M Mset Rset w res r
    FreeEvolvent(0)
    n, result


let rec exec_loop i =
    if i > 5 then 0 
    else
        printfn ""
        printfn "function%O" i
        MySetFunctionNumber i
        let Abounds : List<double> = [GetBoundsA0(); GetBoundsA1()]
        let Bbounds : List<double> = [GetBoundsB0(); GetBoundsB1()]
        let timer = System.Diagnostics.Stopwatch.StartNew()
        let (n: int, res : Result) = GSA_parallel_v2 Abounds Bbounds 2.5 0.0001 10000
        timer.Stop()
        printfn "n = %O" n
        printfn "elapsed = %O" timer.Elapsed
        printfn "F  = %f" res.F
        printfn "Y  = %f %f" (res.Y.Item(0)) (res.Y.Item(1))
        exec_loop (i+1)

let ret = 
    CreateGrishaginFunction()
    exec_loop 1
    FreeGrishaginFunction()

ret