<html>
<head>
<title>Muscle Web Interface</title>
</head>
<body bgcolor="#eef5f5">
<h2>Muscle Web Interface</h2>
<form action='musclewi.py' method='post' enctype="multipart/form-data">
Maximum number of iterations:
<select name="iterat" style="width: 45px" >
  <option value="1" selected="selected">1</option>
  <option value="4">4</option>
  <option value="8">8</option>
  <option value="10">10</option>
  <option value="12">12</option>
  <option value="14">14</option>
  <option value="14">16</option>
</select>
 Output Format:
<select name="output" style="width: 140px" >
  <option value="fasta" selected="selected">FASTA</option>
  <option value="clw">ClustalW2</option>
  <option value="clwstrict">ClustalW2 (Strict)</option>
  <option value="html">HTML</option>
  <option value="msf">MSF</option>
</select>
 Output Order:
 <select name="outorder" style="width: 90px">
  <option value="group" selected="selected">aligned</option>
  <option value="stable">input</option>
 </select>
 <p>Enter or Paste a set of Sequences in any supported format:
 <br/><textarea name="seq" rows="5" cols="90"></textarea><p>
 Or upload a file: <input type="file" name="upfile" />
 <input type='submit' value='Send to Muscle server'></form>
</body></html>
