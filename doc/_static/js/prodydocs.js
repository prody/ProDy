$(document).ready(function() {

  var $window = $(window)
  var items = $("#localtoc a")
  if (items.length) {
    items.prepend('<i class="icon-chevron-right pull-right"></i>')

    var sections = $(".section");
    var section_id = sections[0].id;
    if (section_id.indexOf('.') >= 0) {
      items[0].href = "#" + $(sections[0]).find('>:first-child').attr('id');
    } else {
      items[0].href = "#" + section_id;
    }

    items.first().addClass("navfirst")
    items.last().addClass("navlast")
    $("#localtoc ul").addClass("nav nav-list")
    $("#localtoc > ul").addClass("toc")
  }

  setTimeout(function() {
      $('#sidebar,.subnav').affix({
        offset: {
          top: 140,
          bottom: 20
        }
      });
    }, 100);

  $('.collapse').not('#collapseToolbox').collapse('hide');

  // Alter logo based on URL
  var url = document.URL;
  var logo = document.getElementById("logo");
  if (url.search("plugins") > -1 || url.search("nmwiz") > -1) {
    logo.src = "http://www.csb.pitt.edu/ProDy/_static/nmwiz.png";
  } else if (url.search("evol") > -1 || url.search("panel1-3") > -1) {
    logo.src = "http://www.csb.pitt.edu/ProDy/_static/evol.png";
  }

  var carousel = $('#carousel');
  if (carousel.length) {
    carousel.carousel({'interval': false});
    carousel.on({
        slid: function(){
          var pills = $('.nav-pills > li');
          pills.removeClass('active');
          var index = $('#carousel .item.active').index('#carousel .item');
          var pill = $(pills[index]);
          $(pill).addClass('active');
          if (index == 1) {
            document.getElementById("logo").src = "_static/evol.png";
          } else if (index == 2) {
            document.getElementById("logo").src = "_static/nmwiz.png";
          } else {
            document.getElementById("logo").src = "_static/logo.png";
          }
          window.location.hash = pill.find('a').attr('href').slice(1);
        }});
    // Move carousel based on URL
    if (url.search('evol') > -1) {
      carousel.carousel(1);
    } else if (url.search('nmwiz') > -1) {
      carousel.carousel(2);
    } else if (url.search('downloads') > -1) {
      carousel.carousel(3);
    } else if (url.search('tutorials') > -1) {
      carousel.carousel(4);
    }
  } else {
    // Downloads
    $('tt.download').append(' <i class="icon-download"></i>')
  }

});