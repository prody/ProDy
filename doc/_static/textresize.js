$(document).ready(function() {

  $.cookie.settings = {
      path : "/",
      domain : "",
      expires : 30
  };

  var unit = "pt";
  var size = parseInt($.cookie("ProDyDocFontSize"));
  var defsize = 9;

  if (size) {
      $("body").css("font-size", size + unit);
  } else {
      size = defsize
  }

  function changeFontSize(incr) {
    size = size;
    if (incr) {
      size += incr;
    } else {
      size = defsize;
    }
    $("body").css("font-size", size + unit);
    $.cookie("ProDyDocFontSize", size);
  }

  $('#text-larger').click(function() {changeFontSize(1)});
  $('#text-smaller').click(function() {changeFontSize(-1)});
  $('#text-reset').click(function() {changeFontSize(0)});

});
