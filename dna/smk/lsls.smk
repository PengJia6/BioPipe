rule al:
    input:
         "",
    output:
          "",
    log:
       ""
    benchmark:
             ""
    threads: 8
    run:
        shell("command 1")
        shell("command 2")

rule TODO2:
    input:
         "",
    output:
          "",
    log:
       ""
    benchmark:
             ""
    threads: 8
    shell:
         """
         command 1 
         command 2
         """