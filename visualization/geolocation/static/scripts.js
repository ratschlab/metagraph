var slider = document.getElementById("perc_similarity");
var output = document.getElementById("similarity_cutoff");
output.innerHTML = slider.value; // Display the default slider value

// Update the current slider value (each time you drag the slider handle)
slider.oninput = function() {
    output.innerHTML = this.value;
}
slider.onload = function() {
    output.innerHTML = this.value;
}


function expandTextarea(id) {
    ['input', 'keydown', 'keyup'].forEach(function(e) {
        document.getElementById(id).addEventListener(e, function() {
            this.style.height = 0;
            this.style.height = this.scrollHeight + 'px';
        }, false);
    });
}
expandTextarea('load-from-text');

var num_rows = document.getElementById('matched_samples').rows.length;
var expand_button = document.getElementById('table_toggler');

if (num_rows >= 14) {
    expand_button.style.display = 'inline';

    var show_message = "Show more";
    var hide_message = "Show less";

    expand_button.addEventListener('click', function() {
        if (expand_button.innerHTML == hide_message) {
            expand_button.innerHTML = show_message;
        } else {
            expand_button.innerHTML = hide_message;
        }
    }, false);
}
