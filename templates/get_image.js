var image = document.getElementById("skyplot_img");

document.getElementById("loading_skyplot_img").src = "/static/loading.gif"; 

var downloadingImage = new Image();
downloadingImage.onload = function(){
    image.src = this.src;   
    document.getElementById("loading_skyplot_img").src = "/static/clear.gif"; 
};

var n = current_src.selected.indices[0];

console.log("****", n, current_src.data['cameraratchets'][n])
downloadingImage.src = '/dynamic/image/' + current_src.data['star_ra'][0] + '/' + current_src.data['star_dec'][0] + '/' + current_src.data['cameraratchets'][n] + '/500/300/plot_stars=False';

var rows = document.querySelectorAll('#lc_table_id tr');

// line is zero-based
// line is the row number that you want to see into view after scroll    
var element = rows[n+1];
var topPos = element.offsetTop;

document.getElementById('lc_table_scroll3').scrollTop = topPos - 100;

element.style.backgroundColor = "#AAD";

if (typeof window.lastselectedrow !== 'undefined') {
    if (window.lastselectedrow % 2 == 0) rows[window.lastselectedrow].style.backgroundColor = "#DDD";
    else {rows[window.lastselectedrow].style.backgroundColor = "#FFF";}
}
console.log(window.lastselectedrow, window.lastselectedrow % 2);

window.lastselectedrow = n+1;
