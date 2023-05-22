open System
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
extern void CreateEvovlent()

[<DllImport("evolvent.dll", CallingConvention = CallingConvention.Cdecl)>]
extern void SetBounds(double A0, double A1, double B0, double B1)

[<DllImport("evolvent.dll", CallingConvention = CallingConvention.Cdecl)>]
extern double GetImage0(double x)

[<DllImport("evolvent.dll", CallingConvention = CallingConvention.Cdecl)>]
extern double GetImage1(double x)

[<DllImport("evolvent.dll", CallingConvention = CallingConvention.Cdecl)>]
extern void FreeEvolvent()


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

let minOfPrevAndNewRes (prev : double) (new_p : double) (prev_y : List<double>) (new_y : List<double>) =
    if new_p < prev then (new Result(new_p, new_y))
    else (new Result(prev, prev_y))

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
    if (M = Set.maxElement mset) then 
        let dmp = rmap |> Map.remove (m * (after_new_elem.X - before_new_elem.X) + 
         (after_new_elem.Z - before_new_elem.Z) * (after_new_elem.Z - before_new_elem.Z) / 
         (m * (after_new_elem.X - before_new_elem.X)) - 2.0 * (after_new_elem.Z + before_new_elem.Z))
        let dmp_R = m * (after_new_elem.X - new_elem.X) + (after_new_elem.Z - new_elem.Z) * (after_new_elem.Z - new_elem.Z) / (m * (after_new_elem.X - new_elem.X)) - 2.0 * (after_new_elem.Z + new_elem.Z)
        let dmp_R2 = m * (new_elem.X - before_new_elem.X) + (new_elem.Z - before_new_elem.Z) * (new_elem.Z - before_new_elem.Z) / (m * (new_elem.X - before_new_elem.X)) - 2.0 * (new_elem.Z + before_new_elem.Z)
        let dmp3 = dmp |> Map.add dmp_R after_new_elem
        dmp3 |> Map.add dmp_R2 new_elem, M, m
    else 
        let newM = Set.maxElement mset
        let newm = calc_m newM r
        recalcNewRmap newm w, newM, newm

let mainFunction (PairOfPoints: PairOfPoints) m M mset (rmap: Map<double, Point>) (w : List<Point>) (res : Result) r =
    let x = 0.5 * (PairOfPoints.T.X + PairOfPoints.T_1.X) -  (PairOfPoints.T.Z - PairOfPoints.T_1.Z) / (2.0 * m)
    let new_y = [GetImage0(x);GetImage1(x)]
    let new_point = new Point(x, MyCalculate((new_y.Item(0)), (new_y.Item(1))))
    let new_res : Result = minOfPrevAndNewRes res.F new_point.Z res.Y new_y
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

let rec iter (PairOfPoints: PairOfPoints) epsilon n nmax m M mset rmap w res r = 
    if PairOfPoints.T.X - PairOfPoints.T_1.X <= epsilon or n > nmax then 
        n, res
    else 
        let (PairOfPoints2: PairOfPoints, new_w, newM, new_m, new_rmap, new_mset, new_res) = mainFunction PairOfPoints m M mset rmap w res r
        iter PairOfPoints2 epsilon (n+1) nmax new_m newM new_mset new_rmap new_w new_res r

let GSA (Abounds : List<double>) (Bbounds : List<double>) (r : double) (epsilon : double) (nmax : int) =
    CreateEvovlent()
    SetBounds((Abounds.Item(0)), (Abounds.Item(1)), (Bbounds.Item(0)), (Bbounds.Item(1)))
    let A = 0.0
    let B = 1.0
    let y1 = [GetImage0(A);GetImage1(A)]
    let y2 = [GetImage0(B);GetImage1(B)]
    let point_a = new Point(A, MyCalculate((y1.Item(0)), (y1.Item(1))))
    let point_b = new Point(B, MyCalculate((y2.Item(0)), (y2.Item(1))))
    let w = [point_a; point_b]   
    let M = abs ((point_b.Z - point_a.Z) / (B - A))
    let Mset = Set.empty.Add(M)
    let m = calc_m M r
    let R = m * (B - A) + (point_b.Z - point_a.Z) * (point_b.Z - point_a.Z) / (m * (B - A)) - 2.0 * (point_b.Z + point_a.Z)
    let Rset = Map.empty.Add(R, point_b)
    let res : Result = minOfPrevAndNewRes point_a.Z point_b.Z y1 y2
    iter (new PairOfPoints(point_b, point_a)) epsilon 0 nmax m M Mset Rset w res r

let rec exec_loop i j =
    if i > 5 then 0 
    else
        printfn "function%O" i
        MySetFunctionNumber i
        let Abounds : List<double> = [GetBoundsA0(); GetBoundsA1()]
        let Bbounds : List<double> = [GetBoundsB0(); GetBoundsB1()]
        let timer = System.Diagnostics.Stopwatch.StartNew()
        let (n, res) = GSA Abounds Bbounds 2.5 0.00001 10000
        FreeEvolvent()
        timer.Stop()
        printfn "n = %O" n
        printfn "elapsed = %O" timer.Elapsed
        printfn "F  = %f" res.F
        printfn "Y  = %f %f" (res.Y.Item(0)) (res.Y.Item(1))
        exec_loop (i+1) (j+2)

let ret = 
    CreateGrishaginFunction()
    exec_loop 1 2
    FreeGrishaginFunction()