#assuming lines at infinity or those with x=0 does not occur!
find_gb_lines := proc(f,var::symbol:=x)
  local Lines, ideal:
  global x,y,z,a,b,c,d:
  if var = x then
    Lines := subs([y=a+b*x,z=c+d*x],f):
  elif var = z then #use to get constant x solutions
    Lines := subs([x=c+d*z,y=a+b*z],f):
  fi:
  ideal := [seq(coeff(Lines,var,i),i=0..3)]:
  Groebner:-Basis(ideal, plex(a,b,c,d));  
end proc:

#transforms a cubic form into one in general position (not coaxial with x,y,z and no lines at infinity)
generic_cubic_transform := proc(f,M, homog::boolean:=false)
  local g, G, newV:
  global w,x,y,z:
  
  if homog then g := Groebner:-Homogenize(f,w):
  else g:=f: fi:
  
  newV := M.<w,x,y,z>:
  g := subs([w = newV[1], x=newV[2], y=newV[3], z=newV[4]], g):
  G:=g:
  g := subs(w=1, g):
  g, G, find_gb_lines(g);
end proc:

no_lines := proc(gb):
  factors(gb[1])[2]:
  map(ff->degree(ff[1]),%):
  add(%[i],i=1..nops(%)); #number of lines
end proc:

(* 
Input:  gb, M, eps (optional)
  gb = Groebner basis for computing lines in non-normal form
  M = projective transformation to non-normal form
  eps = positive real number. If given, output is computed numerically with eps specifying the level of precision.
Output: [Lines], [dlines]
  [Lines] = list of lines (original form) parametrized by (s:t)
  [dlines] = list containing 1 or 0 for the lines in Lines, if 1 then lines have multiplicity (singular case) otherwise 0 
*)    
create_lines := proc(gb,M,eps:=0)
  local i,vars, dg, dd, deq,vec, nonzeros, eqns, pairs, pair, mat, facts, line_ds, lines, dlines, Lines, sol, newL, L, s,t:
  global a,b,c,d,w,x,y,z:
  
  vars := [w,x,y,z]:
  
  eqns := gb[2..nops(gb)]:  
  facts := factors(gb[1])[2]:
  line_ds := NULL:
  for deq in facts do:
    if (eps > 0) then line_ds := line_ds, fsolve(deq[1], complex):
    else line_ds := line_ds, solve(deq[1]):
    fi:
  od:
  line_ds := sort([line_ds]): #consistent order of roots

  dg := diff(gb[1],d):
  lines := NULL:
  dlines := NULL: #double or higher "multiplicity" lines
  for dd in line_ds do:
    if (eps>0 and abs(subs(d=dd,dg))<eps) or (eps=0 and evala(subs(d=dd,dg))=0) then
      dlines := dlines, 1:
    else dlines := dlines, 0:
    fi:
    
    solve(subs(d=dd,eqns)):
    lines := lines, subs(%,[w,x,a*w+b*x,c*w+dd*x]):
  od:
  
  Lines := NULL:
  for L in lines do:
    newL := convert(M.<subs([w=s,x=t],L)>, list): #proj. transform back to original
    if (eps=0) then newL:=evala(newL): 
    else newL := evalf(newL): fi:    
    
    #avoiding parameters that are dependent
    nonzeros := NULL:
    for i from 1 to 4 do:
      if newL = 0 then next: fi:
      nonzeros := nonzeros, newL[i] = vars[i]:
    od:
    pairs := combinat:-choose([nonzeros],2);
    for pair in pairs do
      mat := LinearAlgebra:-GenerateMatrix(pair,[s,t]):
      if (eps > 0 and abs(linalg:-det(mat))>eps) or (eps=0 and evala(linalg:-det(mat)) <> 0) then 
        vec := <rhs(pair[1]),rhs(pair[2])>:
        break:
      fi:
    od:
    
    sol := convert(linalg:-inverse(mat).vec,list):
    if (eps = 0) then 
      sol := evala(%): 
      Lines := Lines, evala(subs([s=sol[1],t=sol[2]], newL)):
    else
      Lines := Lines, evalf(subs([s=sol[1],t=sol[2]], newL)):
    fi:    
  od:
  [Lines], [dlines];
end:

(*
Input: Lines
  Lines = is the first output of create_lines
Output: [incs], [poses]
   [incs] = a list of lists of integer indices, such that any element in j in inc[i] satisfy j>i and indicates that Lines[j] intersects with Lines[i]
   [poses] = a list of lists of points in $\mathbb P^3$ (as 4-tuples of algebraic numbers). Such that the k-th point in poses[i] is the point of intersection of Lines[incs[k]] and Lines[i].
*)
compute_incidence := proc(Lines)
  local i, inc, incs, pose, poses, j, sol:
  global w,x,y,z:
  
  incs := NULL:
  poses := NULL:
  for i from 1 to nops(Lines)-1 do:
    inc := NULL:
    pose := NULL:
    for j from i+1 to nops(Lines) do:
      solve(Lines[i]-Lines[j]):
      sol := evala(subs(%,Lines[i])):
      if sol <> [0,0,0,0] then:
        inc := inc, j:
        pose := pose, subs([w=1,x=1,y=1,z=1],sol):
      fi:
    od:
    incs := incs, [inc]:
    poses := poses, [pose]:
  od:

  [incs], [poses]; 
end:

(*
Input: Lines
  Lines = is the first output of create_lines
  eps = positive real number. For numerically computing the incidence relation of the lines, eps specifies the level of precision.
Output: [incs], [poses]
   [incs] = a list of lists of integer indices, such that any element in j in inc[i] satisfy j>i and indicates that Lines[j] intersects with Lines[i]
   [poses] = a list of lists of points in $\mathbb P^3$ (as 4-tuples of algebraic numbers). Such that the k-th point in poses[i] is the point of intersection of Lines[incs[k]] and Lines[i].
*)
compute_incidence_num := proc(Lines, eps)
  local i, inc, incs, pose, poses, j, sol,sub, eqns,vars, L1,L2, err:
  global w,x,y,z:
  
  incs := NULL:
  poses := NULL:
  for i from 1 to nops(Lines)-1 do:  
    inc := NULL:
    pose := NULL:
    for j from i+1 to nops(Lines) do:
      indets(Lines[i]):
      L1 := subs([%[1]=1, %[2]=t1], Lines[i]):
      indets(Lines[j]):
      L2 := subs([%[1]=s2, %[2]=t2], Lines[j]):      
      eqns := L1-L2:
      sub := linalg:-leastsqrs({op(eqns)},{t1,s2,t2}):
      err := linalg:-norm(subs(sub,L1-L2)):
      sol := subs(sub,L1):
      if eps > 0 and err<eps then:
        inc := inc, j:
        pose := pose, sol:
      fi:
    od:
    incs := incs, [inc]:
    poses := poses, [pose]:
  od:

  [incs], [poses]; 
end:

(*
Input: incs, poses. See output of compute_incidence 
Output: [ecks]
  [ecks] = a list of sets, where each set has at least 3 integers. Each integer in such a set are correspond to indices in Lines (see create_lines), indicating which 3 or more lines are concurrent. The common point of intersection of these lines are already encoded in poses.
*)
compute_eckardts := proc(incs, poses)
  local eck, ecks, i, j, k, inc, tested, sol, pose:
  
  ecks := NULL:
  for k from 1 to nops(incs) do:
    inc := incs[k]:
    eck := {}:
    tested := {};  
    for i from 1 to nops(inc)-1 do:  
      if inc[i] in tested then next: fi:
      tested := tested union {inc[i]}:
      eck := {inc[i]}:
      pose := l*~poses[k][i]:
      for j from i+1 to nops(inc) do:
        sol := solve(pose - poses[k][j]):
        if sol <> NULL then #up to constant multiplication
          eck := eck  union {inc[j]}
        fi:    
      od:
      if nops(eck) > 1 then ecks := ecks, (eck union {k}): fi:
    od:
  od:
  [ecks];
end:

