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

(* Input:  gb, M
      gb = Groebner basis for computing lines in non-normal form
      M = projective transformation to non-normal form
   Output: [Lines], [dlines]
    [Lines] = list of lines (original form) parametrized by (s:t)
    [dlines] = list containing 1 or 0 for the lines in Lines, if 1 then lines have multiplicity (singular case) otherwise 0 
*)    
create_lines := proc(gb,M)
  local i,vars, dg, dd, deq,vec, nonzeros, eqns, pairs, pair, mat, facts, line_ds, lines, dlines, Lines, sol, newL, L,s,t:
  global a,b,c,d,w,x,y,z:
  
  vars := [w,x,y,z]:
  
  eqns := gb[2..nops(gb)]:
  facts := factors(gb[1])[2]:
  line_ds := NULL:
  for deq in facts do:
    solve(deq[1]): 
    line_ds := line_ds, %:
  od:
  line_ds := sort([line_ds]): #consistent order of roots

  dg := diff(gb[1],d):
  lines := NULL:
  dlines := NULL: #double or higher "multiplicity" lines
  for dd in line_ds do:
    if evala(subs(d=dd,dg))=0 then
      dlines := dlines, 1:
    else dlines := dlines, 0:
    fi:
    solve(subs(d=dd,eqns)):
    lines := lines, subs(%,[w,x,a*w+b*x,c*w+dd*x]):
  od:
  
  Lines := NULL:
  for L in lines do:
    convert(M.<subs([w=s,x=t],L)>, list): #proj. transform back to original
    newL:=evala(%):  
    
    #avoiding parameters that are dependent
    nonzeros := NULL:
    for i from 1 to 4 do:
      if newL = 0 then next: fi:
      nonzeros := nonzeros, newL[i] = vars[i]:
    od:
    pairs := combinat:-choose([nonzeros],2);
    for pair in pairs do
      mat := LinearAlgebra:-GenerateMatrix(pair,[s,t]):
      if evala(linalg:-det(mat)) <> 0 then 
        vec := <rhs(pair[1]),rhs(pair[2])>:
        break:
      fi:
    od:
    
    convert(linalg:-inverse(mat).vec,list):
    sol := evala(%):
    Lines := Lines, evala(subs([s=sol[1],t=sol[2]], newL)):
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
