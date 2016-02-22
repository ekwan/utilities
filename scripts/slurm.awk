# tells us how busy the cluster is
# usage: squeue2 | awk -f slurm.awk
# pipe the output to sort
# eg: squeue2 | awk -f slurm.awk | sort -n -k3

$2 == "normal" {
    user = $4
    users[user]=1
    if ($5 == "PENDING")
        {
            pending_jobs[user]++
            pending_cores[user] += $10
            total_pending_jobs++
            total_pending_cores += $10
        }
    else if ( $5 == "RUNNING")
        {
            running_jobs[user]++
            running_cores[user] += $10
            total_running_jobs++
            total_running_cores += $10
        }
}

END {
    printf "%d jobs are currently running in the normal partition using %d cores.\n", total_running_jobs, total_running_cores
    printf "%d jobs are pending (%d cores).\n\n", total_pending_jobs, total_pending_cores
    print "user       jobs running      cores in use        jobs pending      cores pending\n"
    for (user in users)
        {
            printf "%15-s  %6d         %6d              %6d           %6d\n", user, running_jobs[user], running_cores[user], pending_jobs[user], pending_cores[user]
        }
}
