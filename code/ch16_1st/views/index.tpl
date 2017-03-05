<!DOCTYPE html>
<html lang="en">
 <head><meta charset="utf-8">
  <title>Vector Filter</title>
  <link href="css/bootstrap.min.css" rel="stylesheet">
 </head>
 <body style="background-color:#e7f5f5;">
 <div class="container"><h2>Vector Filter</h2>
  <form action='/vector_result' method='post'>
   <div class="row">
   <div class="col-sm-8">
   <div class="form-group">
    <label for="seqs">Paste sequence in FASTA format:</label>
    <p><textarea name="seqs" rows="15" cols="80"></textarea></p>
   </div>
   </div>
   </div>
   <div class="row">
   <div class="col-sm-8">
   <div class="form-group">
    <label for="prop">Filter by: </label>
    <div class="radio">
     <label>
      <input type="radio" name="vector" value="customdb">Custom Vectors
     </label>
    </div>
    <div class="radio">
     <label>
      <input type="radio" name="vector" value="ncbivector">NCBI Vector DB
     </label>
    </div>
    <br>
    <button type="submit" class="btn btn-primary">Send
    </button>
   </div>
   </div>
   </div>
  </form>
 </div>
 </body>
</html>
