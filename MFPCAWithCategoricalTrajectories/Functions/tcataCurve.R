tcataCurve=function(df,times=seq(0, 1, by = 0.05),colors=NULL)
{

  # =========================
  # 2. fonction de reconstruction (event -> step function)
  # =========================
  reconstruct_one <- function(df, times){

    df <- df[order(df$time), ]

    etat <- numeric(length(times))
    current <- 0
    j <- 1

    for(i in seq_along(times)){
      while(j <= nrow(df) && df$time[j] <= times[i]){
        current <- df$score[j]
        j <- j + 1
      }
      etat[i] <- current
    }

    data.frame(time = times, score = etat)
  }

  # =========================
  # 3. reconstruction pour tous les sujets / produits / descriptors
  # =========================
  curves <- df %>%
    group_by(descriptor, subject) %>%
    group_split() %>%
    lapply(function(df){

      rec <- reconstruct_one(df, times)

      data.frame(
        descriptor = df$descriptor[1],
        subject = df$subject[1],
        time = rec$time,
        score = rec$score
      )
    }) %>%
    bind_rows()

  # =========================
  # 4. fréquence de citation (moyenne des sujets)
  # =========================
  freq <- curves %>%
    group_by( descriptor, time) %>%
    summarise(freq = mean(score), .groups = "drop")

  # =========================
  # 5. plot final : un panel par produit
  # =========================
  p_curves=ggplot(freq , aes(x = time, y = freq, color = descriptor)) +
    geom_line(linewidth = 1) +
    ylim(0, 1) +
    theme_bw() +
    labs(
      x = "Time",
      y = "Citation frequency",
      color = "Descriptor"
    )+
    scale_color_manual(values = colors)
  return(p_curves)
}