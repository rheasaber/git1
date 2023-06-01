session = null
//QiSession(connected, disconnected, location.host);
QiSession(connected, disconnected, "192.168.1.4");

function connected(s) {
    console.log("Session connected");
    session = s;
    //If you want to subscribe so some events (to send info pepper->tablet) call the function here
}

function disconnected(error) {
    console.log("Session disconnected");
}

function launchEventHome(){
    session.service("ALMemory").then(function (memory) {
            memory.raiseEvent("homeEvent", "paramHome");
        });

   }
function say(sentence){
    session.service("ALAnimatedSpeech").then(function (tts) {
        tts.say(sentence);
        });

   }
function launchEvent2(){
    session.service("ALMemory").then(function (memory) {
            memory.raiseEvent("event2", "param2");
        });

   }
function launchEvent3(){
    session.service("ALMemory").then(function (memory) {
            memory.raiseEvent("event3", "param3");
        });

   }
