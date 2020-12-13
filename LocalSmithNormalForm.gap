#############################################################################
# Smith Normal Form with transformations over the Z/2^N
# written by Alexander Neshitov
#############################################################################

PKG_DIR:="/home/alexander/ch2bt_build/LocalSNF";

SNFTransform := function(mat, N)
  # Smith Normal Form with transformations over Z / 2^N

  # Usage tr:=SNFTransform(mat, N : num_threads = <num_threads>)
  # Arguments:
  #     * mat is matrix over Z,
  #     * N is an integer,
  #     * optional parameter num_thread is the number of parallel threads to be used
  #       in matrix reduction

  # Returns record with fields:
  #   * SNF,
  #   * row_t,
  #   * row_t_inverse,
  #   * col_t,
  #   * rank,
  #   * power
  #
  # That satisfy:
  #  * SNF = row_t * mat * col_t mod 2^N,
  #  * row_t * row_t_inverse = I mod 2^N
  #  * SNF is diagonal, entries are powers of two, every entry divides the next one.
  #  * col_t is invertible over Z/2^N
  #  * rank is the number of nonzero entries of SNF

  local nrows, ncols, num_threads,
        cmd, exec_file,  path, s, mat_str, progress_line,
        tr;

  nrows:= DimensionsMat(mat)[1];
  ncols:= DimensionsMat(mat)[2];
  num_threads := ValueOption("num_threads");
  if num_threads = fail then
    num_threads := 4;
  fi;
  cmd := StringFormatted("src/gap_snf");
  exec_file := Filename(Directory(PKG_DIR), cmd);
  path := DirectoryCurrent();;

  s := InputOutputLocalProcess(path, exec_file, ["-num_threads", String(num_threads)]);
  mat_str := String(mat);;
  RemoveCharacters(mat_str, "[,]");;

  WriteLine(s, Concatenation(StringFormatted("{} {} {}", nrows, ncols, N), mat_str));

  progress_line := ReadLine(s);

  while progress_line<>"result\n" do
    RemoveCharacters(progress_line, "\n");
    Print(Concatenation("\r", progress_line));;
    progress_line := ReadLine(s);

  od;
  tr:=ReadAsFunction(s)();
  CloseStream(s);

  return tr;
end;

SaturationVectors := function(mat, N)
  # Returns columns of row_t_inverse matrix of the SNFTransform
  # corresponding to non-invertible entries of SNF

  local nrows, ncols, num_threads,
        cmd, exec_file, path, s, mat_str, progress_line,
        vectors;

  nrows:= DimensionsMat(mat)[1];
  ncols:= DimensionsMat(mat)[2];
  num_threads := ValueOption("num_threads");
  if num_threads = fail then
    num_threads := 4;
  fi;
  cmd := StringFormatted("src/gap_saturation_vectors");
  exec_file := Filename(Directory(PKG_DIR), cmd);
  path := DirectoryCurrent();;

  s := InputOutputLocalProcess(path, exec_file, ["-num_threads", String(num_threads)]);
  mat_str := String(mat);;
  RemoveCharacters(mat_str, "[,]");;

  WriteLine(s, Concatenation(StringFormatted("{} {} {}", nrows, ncols, N), mat_str));

  progress_line := ReadLine(s);

  while progress_line<>"result\n" do
    RemoveCharacters(progress_line, "\n");
    Print(Concatenation("\r", progress_line));;
    progress_line := ReadLine(s);

  od;
  vectors := ReadAsFunction(s)();
  CloseStream(s);

  return vectors;
end;
