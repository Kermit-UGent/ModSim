Dict(
    "title" => @htl("MODSIM: <strong>Modelling and Simulation of Biosystems I002445</strong>"),

    # # add a disclaimer to the course webpage. Remove it if you dont want to include it.
    "disclaimer" => md"""
    This is the course website for **MODSIM** class tought at Ghent University in the Bachelor of Bioscience Engineering.
    Use it to harness the power of computational thinking.
    """,

    # Highlights the key features of your class to make it more engaging. Remove it if you dont want to include it.
    "highlights" => [
        Dict("name" => "Course Notes", 
             # FIXME Add link to download pdf 
             "text" => md"The MODSIM course notes in pdf format can be found as a download on Ufora or in the (release notes on this GitHub repository)[https://github.com/Kermit-UGent/ModSim/releases].",
             # FIXME change figure to small piece of the cover.
             "img" => "https://github.com/user-attachments/assets/0a178611-f4e3-4572-94a7-ac881aa20f21"
        ),
        Dict("name" => "What is here to be found?!", 
             "text" => md"In this course website you will find info on the course notes, the installation of the software required for this course, and a (static!) view of the exercises and examples. On the left you can find the table of contents with all the pages.",
             # FIXME change picture
             "img" => "https://github.com/user-attachments/assets/796b2e07-b741-49cd-8ae9-5ce8d9fff4d2"
        ),
        Dict("name" => "Questions, office hours, errata and bugs",
             # FIXME check content 
             "text" => md"For question you can check the fora on Ufora, ask the teaching responsible or assistents in class. Further, there are office hours on XYZ day from AAA to BBB in office number XYZ.123. Errata will be anounced on the errate page If you find something that looks wonky, is misspelledt (pun intended) or is bugged: please open a bug report here (https://github.com/Kermit-UGent/ModSim/issues/) or report it to one of the teaching assistents.",
             # FIXME change picture
             "img" => "https://github.com/user-attachments/assets/56e59af1-8566-40ec-85ac-7184926a29fc"
        ),
        Dict("name" => "Interactivity with Julia and Pluto.jl",
             "text" => md"""
             Thanks to Julia and Pluto.jl, this website and the exercises accompanying the course notes are built using real code, and instead of a book, we have a series of interactive notebooks.
             **On the website, students can play with sliders, buttons and images to interact with our simulations.**
             You can even go further, and modify and run any code on our website!
             """,
             "img" => "https://user-images.githubusercontent.com/6933510/136196607-16207911-53be-4abb-b90e-d46c946e6aaf.gif"
             ),
     #    Dict("name" => "Learn Julia",
     #         "text" => md"""
     #         In literature it's not enough to just know the technicalities of grammar.
     #         In music it's not enough to learn the scales. The goal is to communicate experiences and emotions.
     #         For a computer scientist, it's not enough to write a working program,
     #         the program should be written with beautiful high level abstractions that speak to your audience.
     #         **Julia is designed with this purpose in mind, use it in your teaching to harness its power.**
     #         """,
     #         "img" => "https://user-images.githubusercontent.com/6933510/136203632-29ce0a96-5a34-46ad-a996-de55b3bcd380.png"
     #    )
    ]
)