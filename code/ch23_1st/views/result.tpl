<html lang="en">
  <head>
    <meta charset="utf-8">
    <title>Muscle Web Interface</title>
    <link href="css/bootstrap.min.css" rel="stylesheet">
  </head>
  <body>
  <div class="container">
      % if bad_option:
      <h1>Bad Request</h1>
      Use the options provided in the form. Error in {{bad_option}}
      % elif not bad_option:
      {{!result_output}}
      % end
      <p>Go back to the <a href="/">home page</a></p>
  </div>
  </body>
</html>
