
<!-- This is the project specific website template -->
<!-- It can be changed as liked or replaced by other content -->

<?php

$domain=ereg_replace('[^\.]*\.(.*)$','\1',$_SERVER['HTTP_HOST']);
$group_name=ereg_replace('([^\.]*)\..*$','\1',$_SERVER['HTTP_HOST']);
$themeroot='r-forge.r-project.org/themes/rforge/';

echo '<?xml version="1.0" encoding="UTF-8"?>';
?>
<!DOCTYPE html
	PUBLIC "-//W3C//DTD XHTML 1.0 Transitional//EN"
	"http://www.w3.org/TR/xhtml1/DTD/xhtml1-transitional.dtd">
<html xmlns="http://www.w3.org/1999/xhtml" xml:lang="en" lang="en   ">

  <head>
	<meta http-equiv="Content-Type" content="text/html; charset=UTF-8" />
	<title><?php echo $group_name; ?></title>
	<link href="http://<?php echo $themeroot; ?>styles/estilo1.css" rel="stylesheet" type="text/css" />
  </head>

<body>

<!-- R-Forge Logo -->
<table border="0" width="100%" cellspacing="0" cellpadding="0">
<tr><td>
<a href="http://r-forge.r-project.org/"><img src="http://<?php echo $themeroot; ?>/imagesrf/logo.png" border="0" alt="R-Forge Logo" /> </a> </td> </tr>
</table>


<!-- get project title  -->
<!-- own website starts here, the following may be changed as you like -->

<?php if ($handle=fopen('http://'.$domain.'/export/projtitl.php?group_name='.$group_name,'r')){
$contents = '';
while (!feof($handle)) {
	$contents .= fread($handle, 8192);
}
fclose($handle);
echo $contents; } ?>

<!-- end of project description -->

<hr />

<h2>Description</h2>

<p> The <a href="http://<?php echo $domain; ?>/projects/<?php echo $group_name; ?>/">R-Forge project <tt>ebuzh</tt></a> currently hosts a single R package: <strong>biostatUZH</strong>. </p>

<p> For information about the package contents see the
<a href="https://<?php echo $domain; ?>/scm/viewvc.php/*checkout*/pkg/biostatUZH/DESCRIPTION?root=<?php echo $group_name; ?>&content-type=text%2Fplain">DESCRIPTION file</a>
and the
<a href="https://<?php echo $domain; ?>/scm/viewvc.php/*checkout*/pkg/biostatUZH/inst/NEWS.Rd?root=<?php echo $group_name; ?>&content-type=text%2Fplain">latest NEWS</a>.</p>

<h2>Installation</h2>

<p> You can easily install the current package version in <tt>R</tt>. </p>

<h3>Windows &amp; Linux</h3>

<div style="text-align:center; font-size:x-large; font-family:monospace">
install.packages("biostatUZH", repos = "http://R-Forge.R-project.org")
</div>

<h3>Mac OS X</h3>

<p> Since R-Forge does no longer offer binaries for MacOS X,
it is necessary to specify the additional argument <tt>type = "source"</tt>: </p>
<div style="text-align:center; font-size:x-large; font-family:monospace">
install.packages("biostatUZH", repos = "http://R-Forge.R-project.org", type = "source")
</div>

</body>
</html>
